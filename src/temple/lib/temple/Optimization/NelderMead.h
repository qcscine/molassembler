/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#ifndef INCLUDE_TEMPLE_OPTIMIZATION_NELDER_MEAD_H
#define INCLUDE_TEMPLE_OPTIMIZATION_NELDER_MEAD_H

#include <Eigen/Core>
#include "temple/Functional.h"
#include "temple/TinySet.h"

namespace temple {

template<typename FloatType = double>
struct NelderMead {
  using MatrixType = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;

  //! Type returned from an optimization
  struct OptimizationReturnType {
    //! Number of iterations
    unsigned iterations;
    //! Final function value
    FloatType value;
  };

  template<
    typename UpdateFunction,
    typename Checker
  > static OptimizationReturnType minimize(
    Eigen::Ref<MatrixType> vertices,
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

    const unsigned N = vertices.rows();
    assert(vertices.cols() == N + 1);

    struct IndexValuePair {
      unsigned column;
      FloatType value;

      bool operator < (const IndexValuePair& other) const {
        return value < other.value;
      }
    };

    std::vector<IndexValuePair> values = temple::sort(
      temple::map(
        temple::iota<unsigned>(N + 1),
        [&](const unsigned i) -> IndexValuePair {
          return {
            i,
            function(vertices.col(i))
          };
        }
      )
    );

    auto replaceWorst = [&values](
      const VectorType& vertex,
      const FloatType value,
      Eigen::Ref<MatrixType> verticesRef
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
      verticesRef.col(replacementPair.column) = vertex;
    };

    FloatType standardDeviation;
    unsigned iteration = 0;
    do {
      const VectorType simplexCentroid = vertices.rowwise().sum() / (N + FloatType {1});

      /* Reflect */
      const VectorType reflectedVertex = (
        simplexCentroid
        + alpha * (simplexCentroid - vertices.col(values.back().column))
      );
      const FloatType reflectedValue = function(reflectedVertex);

      if(
        values.front().value <= reflectedValue
        && reflectedValue < values.at(N - 1).value
      ) {
        replaceWorst(reflectedVertex, reflectedValue, vertices);
      } else if(reflectedValue < values.front().value) {
        /* Expansion */
        const VectorType expandedVertex = simplexCentroid + gamma * (reflectedVertex - simplexCentroid);
        const FloatType expandedValue = function(expandedVertex);

        if(expandedValue < reflectedValue) {
          // Replace the worst value with the expanded point
          replaceWorst(expandedVertex, expandedValue, vertices);
        } else {
          // Replace the worst value with the reflected point
          replaceWorst(reflectedVertex, reflectedValue, vertices);
        }
      } else {
        /* Contraction */
        const VectorType contractedVertex = simplexCentroid + rho * (vertices.col(values.back().column) - simplexCentroid);
        const FloatType contractedValue = function(contractedVertex);
        if(contractedValue < values.back().value) {
          replaceWorst(contractedVertex, contractedValue, vertices);
        } else {
          /* Shrink */
          const unsigned bestColumn = values.front().column;
          const auto& bestVertex = vertices.col(bestColumn);
          for(unsigned i = 0; i < bestColumn; ++i) {
            vertices.col(i) = bestVertex + sigma * (vertices.col(i) - bestVertex);
          }
          for(unsigned i = bestColumn + 1; i <= N; ++i) {
            vertices.col(i) = bestVertex + sigma * (vertices.col(i) - bestVertex);
          }
        }
      }

      // Calculate standard deviation of values
      FloatType sum = 0;
      for(const auto& p : values) {
        sum += p.value;
      }
      const FloatType average = sum / (N + 1);
      standardDeviation = 0;
      for(const auto& p : values) {
        const FloatType diff = p.value - average;
        standardDeviation += diff * diff;
      }
      standardDeviation = std::sqrt(standardDeviation / (N + 1));
      ++iteration;
    } while(check.shouldContinue(iteration, standardDeviation));

    return {
      iteration,
      values.front().value
    };
  }
};

} // namespace temple

#endif
