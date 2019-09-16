/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Nelder-Mead over the SO(3) manifold
 */

#ifndef INCLUDE_TEMPLE_OPTIMIZATION_SO3_NELDER_MEAD_H
#define INCLUDE_TEMPLE_OPTIMIZATION_SO3_NELDER_MEAD_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include "temple/Adaptors/AllPairs.h"
#include "temple/Functional.h"
#include "temple/TinySet.h"

namespace temple {

template<typename FloatType = double>
struct SO3NelderMead {
  using Matrix = Eigen::Matrix<FloatType, 3, 3>;
  using Parameters = std::array<Matrix, 4>;

  //! Type returned from an optimization
  struct OptimizationReturnType {
    //! Number of iterations
    unsigned iterations;
    //! Final function value
    FloatType value;
  };

  static FloatType distanceSquared(const Matrix& a, const Matrix& b) {
    return FloatType {0.5} * (a.transpose() * b).log().squaredNorm();
  }

  static Matrix geodesicExtrapolation(const Matrix& a, const Matrix& b, const FloatType tau) {
    return b * (tau * (b.transpose() * a).log()).exp();
  }

  static bool isOrthogonal(const Matrix& m) {
    return (m * m.transpose()).isApprox(Matrix::Identity(), 1e-5);
  }

  static Matrix karcherMean(const Parameters& points) {
    auto calculateOmega = [](const Parameters& p, const Matrix& speculativeMean) {
      Matrix omega = Matrix::Zero();
      for(const Matrix& point : p) {
        Matrix U = (speculativeMean.transpose() * point).log();
        omega += 0.5 * (U - U.transpose());
      }
      omega /= 4;
      return omega;
    };

    constexpr FloatType delta = 1e-5;
    Matrix q = points.front();
    Matrix omega = calculateOmega(points, q);

    unsigned iterations = 0;
    while(omega.norm() >= delta && iterations < 100) {
      q = q * omega.exp();
      assert(q.allFinite());
      assert(isOrthogonal(q));
      omega = calculateOmega(points, q);
      ++iterations;
    }

    return q;
  }

  static Matrix randomRotation() {
    auto A = Matrix::Random();
    Eigen::ColPivHouseholderQR<Matrix> decomposition(A);
    Matrix Q = decomposition.householderQ();
    auto R = decomposition.matrixR();

    Matrix intermediate = Matrix::Zero();
    for(unsigned i = 0; i < 3; ++i) {
      double value = R(i, i);
      if(value < 0) {
        intermediate(i, i) = -1;
      } else if(value > 0) {
        intermediate(i, i) = 1;
      }
    }
    Q = Q * intermediate;

    // Now Q is in O(n), but not yet in SO(n), which we can ensure with:
    if(Q.determinant() < 0) {
      Q.col(0).swap(Q.col(1));
    }

    // Is it really orthogonal?
    assert(isOrthogonal(Q));
    assert(Q.allFinite());
    return Q;
  }

  static Parameters randomParameters() {
    constexpr FloatType ballRadiusSquared = M_PI * M_PI;
    Parameters parameters;
    parameters.front() = randomRotation();
    for(unsigned i = 1; i < 4; ++i) {
      Matrix R;
      do {
        R = randomRotation();
      } while(
        temple::any_of(
          temple::iota<unsigned>(i),
          [&](const unsigned j) -> bool {
            return distanceSquared(R, parameters.at(j)) >= ballRadiusSquared;
          }
        )
      );
      parameters.at(i) = R;
    }

    return parameters;
  }

  template<
    typename UpdateFunction,
    typename Checker
  > static OptimizationReturnType minimize(
    Parameters& points,
    UpdateFunction&& function,
    Checker&& check
  ) {
    constexpr FloatType alpha = 1; // Reflection coefficient
    constexpr FloatType gamma = 2; // Expansion coefficient
    constexpr FloatType rho = 0.5; // Contraction coefficient
    constexpr FloatType sigma = 0.5; // Shrink coefficient

    static_assert(0 < alpha, "Alpha bounds not met");
    static_assert(1 < gamma, "Gamma bounds not met");
    static_assert(0 < rho && rho <= 0.5, "Rho bounds not met");

    struct IndexValuePair {
      unsigned column;
      FloatType value;

      bool operator < (const IndexValuePair& other) const {
        return value < other.value;
      }
    };

    /* We need the points of the simplex to always lie within a ball of radius
     * pi/2 so that geodesics and the karcher mean are unique.
     *
     * We think that if all pairs of points have distance less than pi from one
     * another, they are within a ball of radius pi/2.
     */
    constexpr FloatType ballRadiusSquared = M_PI * M_PI;
    if(
      temple::any_of(
        temple::adaptors::allPairs(points),
        [](const Matrix& p1, const Matrix& p2) -> bool {
          return distanceSquared(p1, p2) >= ballRadiusSquared;
        }
      )
    ) {
      throw std::logic_error(
        "Initial simplex points do not lie within ball of radius pi/2"
      );
    }

    // Sort the vertex values
    std::vector<IndexValuePair> values = temple::sort(
      temple::map(
        temple::iota<unsigned>(4),
        [&](const unsigned i) -> IndexValuePair {
          return {
            i,
            function(points.at(i))
          };
        }
      )
    );

    auto replaceWorst = [&values](
      const Matrix& point,
      const FloatType value,
      Parameters& simplexVertices
    ) {
      IndexValuePair replacementPair {
        values.back().column,
        value
      };
      // Insert the replacement into values, keeping ordering
      temple::TinySet<IndexValuePair>::checked_insert(values, replacementPair);
      // Drop the worst value
      values.pop_back();
      // Replace the worst column
      simplexVertices.at(replacementPair.column) = point;
    };

    auto ballCheckingFunction = [](
      auto&& objectiveFunction,
      const Parameters& simplexVertices,
      const Matrix& speculativePoint,
      const unsigned replacingIndex
    ) -> FloatType {
      for(unsigned i = 0; i < 4; ++i) {
        if(i == replacingIndex) {
          continue;
        }

        if(distanceSquared(speculativePoint, simplexVertices.at(i)) >= ballRadiusSquared) {
          /* The new point will lie outside a ball of radius pi / 2 for the
           * existing points. Geodesics may no longer be unique and the Karcher
           * mean may be incalculable. Returning a near-infinite function value
           * will discourage the use of this point.
           */
          return std::numeric_limits<FloatType>::max();
        }
      }

      return objectiveFunction(speculativePoint);
    };

    FloatType standardDeviation;
    unsigned iteration = 0;
    do {
      const Matrix simplexCentroid = karcherMean(points);
      assert(isOrthogonal(simplexCentroid));

      /* Reflect */
      const Matrix reflectedVertex = geodesicExtrapolation(points.at(values.back().column), simplexCentroid, -alpha);
      const FloatType reflectedValue = ballCheckingFunction(function, points, reflectedVertex, values.back().column);

      if(
        values.front().value <= reflectedValue
        && reflectedValue < values.at(2).value
      ) {
        replaceWorst(reflectedVertex, reflectedValue, points);
      } else if(reflectedValue < values.front().value) {
        /* Expansion */
        const Matrix expandedVertex = geodesicExtrapolation(points.at(values.back().column), simplexCentroid, -gamma);
        const FloatType expandedValue = ballCheckingFunction(function, points, expandedVertex, values.back().column);

        if(expandedValue < reflectedValue) {
          // Replace the worst value with the expanded point
          replaceWorst(expandedVertex, expandedValue, points);
        } else {
          // Replace the worst value with the reflected point
          replaceWorst(reflectedVertex, reflectedValue, points);
        }
      } else {
        /* Contraction */
        const Matrix contractedVertex = geodesicExtrapolation(points.at(values.back().column), simplexCentroid, rho);
        const FloatType contractedValue = function(contractedVertex);
        // NOTE: No need to worry about ball radius in inside contraction
        if(contractedValue < values.back().value) {
          replaceWorst(contractedVertex, contractedValue, points);
        } else {
          /* Shrink */
          const Matrix& bestVertex = points.at(values.front().column);

          // Shrink all points besides the best one and recalculate function values
          for(unsigned i = 1; i < 4; ++i) {
            auto& value = values.at(i);
            auto& vertex = points.at(value.column);
            vertex = geodesicExtrapolation(vertex, bestVertex, sigma);
            value.value = function(vertex);
            // NOTE: No need to worry about ball radius in shrink operation
          }
          temple::inplace::sort(values);
        }
      }

      // Calculate standard deviation of values
      const FloatType average = temple::accumulate(
        values,
        FloatType {0},
        [](const FloatType carry, const IndexValuePair& pair) -> FloatType {
          return carry + pair.value;
        }
      ) / 4;
      standardDeviation = std::sqrt(
        temple::accumulate(
          values,
          FloatType {0},
          [average](const FloatType carry, const IndexValuePair& pair) -> FloatType {
            const FloatType diff = pair.value - average;
            return carry + diff * diff;
          }
        ) / 4
      );
      ++iteration;
    } while(check.shouldContinue(iteration, values.front().value, standardDeviation));

    return {
      iteration,
      values.front().value
    };
  }
};

} // namespace temple

#endif
