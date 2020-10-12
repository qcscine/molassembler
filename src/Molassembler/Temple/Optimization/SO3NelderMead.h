/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Nelder-Mead over the SO(3) manifold
 */

#ifndef INCLUDE_TEMPLE_OPTIMIZATION_SO3_NELDER_MEAD_H
#define INCLUDE_TEMPLE_OPTIMIZATION_SO3_NELDER_MEAD_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Functional.h"

namespace Scine {
namespace Molassembler {
namespace Temple {

/**
 * @brief Nelder-Mead optimization on SO(3) manifold
 *
 * @tparam FloatType Type to represent floating point numbers
 */
template<typename FloatType = double>
struct SO3NelderMead {
  using Matrix = Eigen::Matrix<FloatType, 3, 3>;
  //! Four 3x3 matrices form the simplex vertices
  struct Parameters {
    Eigen::Matrix<FloatType, 3, 12> matrix;

    decltype(auto) at(const unsigned i) {
      assert(i < 4);
      return matrix.template block<3, 3>(0, 3 * i);
    }

    decltype(auto) at(const unsigned i) const {
      assert(i < 4);
      return matrix.template block<3, 3>(0, 3 * i);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

  //! Type returned from an optimization
  struct OptimizationReturnType {
    //! Number of iterations
    unsigned iterations;
    //! Final function value
    FloatType value;
    //! Simplex vertex with minimal function value
    unsigned minimalIndex;
  };

  struct Manifold {
    // Returns the skew-symmetric parts of m
    template<typename Derived>
    static Matrix skew(const Eigen::MatrixBase<Derived>& m) {
      return 0.5 * (m - m.transpose());
    }

    template<typename DerivedA, typename DerivedB>
    static Matrix log(const Eigen::MatrixBase<DerivedA>& X, const Eigen::MatrixBase<DerivedB>& Y) {
      return skew((X.transpose() * Y).log().real());
    }

    template<typename DerivedA, typename DerivedB>
    static Matrix exp(const Eigen::MatrixBase<DerivedA>& X, const Eigen::MatrixBase<DerivedB>& Y) {
      return X * (Y.exp());
    }

    template<typename DerivedA, typename DerivedB>
    static FloatType distanceSquared(const Eigen::MatrixBase<DerivedA>& X, const Eigen::MatrixBase<DerivedB>& Y) {
      return FloatType {0.5} * log(X, Y).squaredNorm();
    }

    template<typename DerivedA, typename DerivedB>
    static Matrix geodesic(const Eigen::MatrixBase<DerivedA>& a, const Eigen::MatrixBase<DerivedB>& b, const FloatType tau) {
      return exp(b, tau * log(b, a));
    }

    template<typename Derived>
    static bool contains(const Eigen::MatrixBase<Derived>& m) {
      return (m * m.transpose()).isApprox(Matrix::Identity(), 1e-5);
    }

    static Matrix karcherMean(const Parameters& points, const unsigned excludeIdx) {
      auto calculateOmega = [excludeIdx](const Parameters& p, const Matrix& speculativeMean) {
        Matrix omega = Matrix::Zero();
        for(unsigned i = 0; i < 4; ++i) {
          if(i == excludeIdx) {
            continue;
          }

          omega += Manifold::log(speculativeMean, p.at(i));
        }
        omega /= 4;
        return omega;
      };

      constexpr FloatType delta = 1e-5;
      Matrix q = (excludeIdx == 0) ? points.at(1) : points.at(0);
      Matrix omega = calculateOmega(points, q);

      unsigned iterations = 0;
      while(omega.norm() >= delta && iterations < 100) {
        q = Manifold::exp(q, omega);
        assert(q.allFinite());
        assert(contains(q));
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
      assert(Q.allFinite());
      assert(contains(Q));
      return Q;
    }
  };

  static Parameters randomParameters() {
    constexpr FloatType ballRadiusSquared = M_PI * M_PI;
    Parameters parameters;
    parameters.at(0) = Manifold::randomRotation();
    for(unsigned i = 1; i < 4; ++i) {
      Matrix R;
      do {
        R = Manifold::randomRotation();
      } while(
        Temple::any_of(
          Temple::iota<unsigned>(i),
          [&](const unsigned j) -> bool {
            return Manifold::distanceSquared(R, parameters.at(j)) >= ballRadiusSquared;
          }
        )
      );
      parameters.at(i) = R;
    }

    return parameters;
  }

  struct IndexValuePair {
    unsigned column;
    FloatType value;

    bool operator < (const IndexValuePair& other) const {
      return value < other.value;
    }
  };

  static FloatType valueStandardDeviation(const std::vector<IndexValuePair>& sortedPairs) {
    const unsigned V = sortedPairs.size();
    // Calculate standard deviation of values
    const FloatType average = Temple::accumulate(
      sortedPairs,
      FloatType {0},
      [](const FloatType carry, const IndexValuePair& pair) -> FloatType {
        return carry + pair.value;
      }
    ) / V;
    return std::sqrt(
      Temple::accumulate(
        sortedPairs,
        FloatType {0},
        [average](const FloatType carry, const IndexValuePair& pair) -> FloatType {
          const FloatType diff = pair.value - average;
          return carry + diff * diff;
        }
      ) / V
    );
  }

  template<typename UpdateFunction>
  static void shrink(
    Parameters& points,
    std::vector<IndexValuePair>& values,
    UpdateFunction&& function
  ) {
    constexpr FloatType shrinkCoefficient = 0.5;
    static_assert(0 < shrinkCoefficient && shrinkCoefficient < 1, "Shrink coefficient bounds not met");
    const Matrix& bestVertex = points.at(values.front().column);

    // Shrink all points besides the best one and recalculate function values
    for(unsigned i = 1; i < 4; ++i) {
      auto& value = values.at(i);
      Eigen::Ref<Eigen::Matrix3d> vertex = points.at(value.column);
      vertex = Manifold::geodesic(vertex, bestVertex, shrinkCoefficient);
      value.value = function(vertex);
      // NOTE: No need to worry about ball radius in shrink operation
    }
    Temple::sort(values);
  }

  static void replaceWorst(
    std::vector<IndexValuePair>& sortedPairs,
    const Matrix& newVertex,
    const FloatType newValue,
    Parameters& vertices
  ) {
    assert(std::is_sorted(std::begin(sortedPairs), std::end(sortedPairs)));
    IndexValuePair replacementPair {
      sortedPairs.back().column,
      newValue
    };
    sortedPairs.insert(
      std::lower_bound(
        std::begin(sortedPairs),
        std::end(sortedPairs),
        replacementPair
      ),
      replacementPair
    );
    // Drop the worst value
    sortedPairs.pop_back();
    // Replace the worst column
    vertices.at(replacementPair.column) = newVertex;
  }

  template<
    typename UpdateFunction,
    typename Checker
  > static OptimizationReturnType minimize(
    Parameters& points,
    UpdateFunction&& function,
    Checker&& check
  ) {
    constexpr FloatType reflectionCoefficient = 1;
    constexpr FloatType expansionCoefficient = 2;
    constexpr FloatType contractionCoefficient = 0.5;

    static_assert(0 < reflectionCoefficient, "Reflection coefficient bounds not met");
    static_assert(1 < expansionCoefficient, "Expansion coefficient bounds not met");
    static_assert(0 < contractionCoefficient && contractionCoefficient <= 0.5, "Contraction coefficient bounds not met");

    /* We need the points of the simplex to always lie within a ball of radius
     * pi/2 so that geodesics and the karcher mean are unique.
     *
     * We think that if all pairs of points have distance less than pi from one
     * another, they are within a ball of radius pi/2.
     */
    constexpr FloatType ballRadiusSquared = M_PI * M_PI;
    if(
      Temple::any_of(
        Temple::Adaptors::allPairs(Temple::iota<unsigned>(4)),
        [&points](const unsigned i, const unsigned j) -> bool {
          return Manifold::distanceSquared(points.at(i), points.at(j)) >= ballRadiusSquared;
        }
      )
    ) {
      throw std::logic_error(
        "Initial simplex points do not lie within ball of radius pi/2"
      );
    }

    // Sort the vertex values
    std::vector<IndexValuePair> values = Temple::sorted(
      Temple::map(
        Temple::iota<unsigned>(4),
        [&](const unsigned i) -> IndexValuePair {
          return {
            i,
            function(points.at(i))
          };
        }
      )
    );
    assert(values.size() == 4);

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

        if(Manifold::distanceSquared(speculativePoint, simplexVertices.at(i)) >= ballRadiusSquared) {
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
      const Matrix simplexCentroid = Manifold::karcherMean(points, values.back().column);
      const Matrix& worstVertex = points.at(values.back().column);
      const FloatType worstVertexValue = values.back().value;
      const FloatType bestVertexValue = values.front().value;

      /* Reflect */
      const Matrix reflectedVertex = Manifold::geodesic(worstVertex, simplexCentroid, -reflectionCoefficient);
      const FloatType reflectedValue = ballCheckingFunction(function, points, reflectedVertex, values.back().column);

      if(reflectedValue < bestVertexValue) {
        /* Expansion */
        const Matrix expandedVertex = Manifold::geodesic(worstVertex, simplexCentroid, -expansionCoefficient);
        const FloatType expandedValue = ballCheckingFunction(function, points, expandedVertex, values.back().column);

        if(expandedValue < reflectedValue) {
          // Replace the worst value with the expanded point
          replaceWorst(values, expandedVertex, expandedValue, points);
        } else {
          // Replace the worst value with the reflected point
          replaceWorst(values, reflectedVertex, reflectedValue, points);
        }
      } else if(bestVertexValue <= reflectedValue && reflectedValue < values.at(2).value) {
        replaceWorst(values, reflectedVertex, reflectedValue, points);
      } else if(values.at(2).value <= reflectedValue && reflectedValue < worstVertexValue) {
        /* Outside contraction */
        const Matrix outsideContractedVertex = Manifold::geodesic(worstVertex, simplexCentroid, -(reflectionCoefficient * contractionCoefficient));
        const FloatType outsideContractedValue = ballCheckingFunction(function, points, reflectedVertex, values.back().column);
        if(outsideContractedValue <= reflectedValue) {
          replaceWorst(values, outsideContractedVertex, outsideContractedValue, points);
        } else {
          shrink(points, values, function);
        }
      } else {
        /* Inside contraction */
        const Matrix insideContractedVertex = Manifold::geodesic(worstVertex, simplexCentroid, contractionCoefficient);
        const FloatType insideContractedValue = function(insideContractedVertex);
        // NOTE: No need to worry about ball radius in inside contraction
        if(insideContractedValue < worstVertexValue) {
          replaceWorst(values, insideContractedVertex, insideContractedValue, points);
        } else {
          shrink(points, values, function);
        }
      }

      standardDeviation = valueStandardDeviation(values);
      ++iteration;
    } while(check.shouldContinue(iteration, values.front().value, standardDeviation));

    return {
      iteration,
      values.front().value,
      values.front().column
    };
  }
};

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
