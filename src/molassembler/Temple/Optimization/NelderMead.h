/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Nelder-Mead in Euclidean space
 */

#ifndef INCLUDE_TEMPLE_OPTIMIZATION_NELDER_MEAD_H
#define INCLUDE_TEMPLE_OPTIMIZATION_NELDER_MEAD_H

#include <Eigen/Core>
#include "temple/Functional.h"

namespace Scine {
namespace temple {

/**
 * @brief Nelder-Mead optimization
 *
 * @tparam FloatType Type to represent floating point numbers
 */
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
    //! Simplex vertex with minimal function value
    unsigned minimalIndex;
  };

  static VectorType generateVertex(
    const FloatType coefficient,
    const VectorType& centroid,
    const VectorType& worstVertex
  ) {
    return centroid + coefficient * (centroid - worstVertex);
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
    const FloatType average = temple::accumulate(
      sortedPairs,
      FloatType {0},
      [](const FloatType carry, const IndexValuePair& pair) -> FloatType {
        return carry + pair.value;
      }
    ) / V;
    return std::sqrt(
      temple::accumulate(
        sortedPairs,
        FloatType {0},
        [average](const FloatType carry, const IndexValuePair& pair) -> FloatType {
          const FloatType diff = pair.value - average;
          return carry + diff * diff;
        }
      ) / V
    );
  }

  //! Calculates the simplex centroid excluding the worst vertex
  static VectorType centroid(
    const MatrixType& vertices,
    const std::vector<IndexValuePair>& pairs
  ) {
    return (vertices.rowwise().sum() - vertices.col(pairs.back().column)) / (vertices.cols() - 1);
  }

  template<typename UpdateFunction>
  static void shrink(
    Eigen::Ref<MatrixType> vertices,
    std::vector<IndexValuePair>& pairs,
    UpdateFunction&& function
  ) {
    constexpr FloatType shrinkCoefficient = 0.5;
    static_assert(0 < shrinkCoefficient && shrinkCoefficient < 1, "Shrink coefficient bounds not met");

    const unsigned V = vertices.cols();
    const unsigned bestColumn = pairs.front().column;
    const auto& bestVertex = vertices.col(bestColumn);

    // Shrink all points besides the best one and recalculate function values
    for(unsigned i = 1; i < V; ++i) {
      auto& value = pairs.at(i);
      auto vertex = vertices.col(value.column);
      vertex = bestVertex + shrinkCoefficient * (vertex - bestVertex);
      value.value = function(vertex);
    }
    temple::inplace::sort(pairs);
  }

  static void replaceWorst(
    std::vector<IndexValuePair>& sortedPairs,
    const VectorType& newVertex,
    const FloatType newValue,
    Eigen::Ref<MatrixType> vertices
  ) {
    assert(std::is_sorted(std::begin(sortedPairs), std::end(sortedPairs)));
    IndexValuePair replacementPair {
      sortedPairs.back().column,
      newValue
    };
    // Insert the replacement into values, keeping ordering
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
    vertices.col(replacementPair.column) = newVertex;
  }

  template<
    typename UpdateFunction,
    typename Checker
  > static OptimizationReturnType minimize(
    Eigen::Ref<MatrixType> vertices,
    UpdateFunction&& function,
    Checker&& check
  ) {
    constexpr FloatType reflectionCoefficient = 1;
    constexpr FloatType expansionCoefficient = 2;
    constexpr FloatType contractionCoefficient = 0.5;

    static_assert(0 < reflectionCoefficient, "Reflection coefficient bounds not met");
    static_assert(1 < expansionCoefficient, "Expansion coefficient bounds not met");
    static_assert(0 < contractionCoefficient && contractionCoefficient <= 0.5, "Contraction coefficient bounds not met");

    const unsigned N = vertices.rows();
    assert(vertices.cols() == N + 1);

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

    FloatType standardDeviation;
    unsigned iteration = 0;
    do {
      const VectorType simplexCentroid = centroid(vertices, values);
      const VectorType& worstVertex = vertices.col(values.back().column);
      const FloatType worstVertexValue = values.back().value;
      const FloatType bestVertexValue = values.front().value;

      /* Reflect */
      const VectorType reflectedVertex = generateVertex(reflectionCoefficient, simplexCentroid, worstVertex);
      const FloatType reflectedValue = function(reflectedVertex);

      if(reflectedValue < bestVertexValue) {
        /* Expansion */
        const VectorType expandedVertex = generateVertex(expansionCoefficient, simplexCentroid, worstVertex);
        const FloatType expandedValue = function(expandedVertex);

        if(expandedValue < reflectedValue) {
          // Replace the worst value with the expanded point
          replaceWorst(values, expandedVertex, expandedValue, vertices);
        } else {
          // Replace the worst value with the reflected point
          replaceWorst(values, reflectedVertex, reflectedValue, vertices);
        }
      } else if(bestVertexValue <= reflectedValue && reflectedValue < values.at(N - 1).value) {
        replaceWorst(values, reflectedVertex, reflectedValue, vertices);
      } else if(values.at(N - 1).value <= reflectedValue && reflectedValue < worstVertexValue) {
        /* Outside contraction */
        const VectorType outsideContractedVertex = generateVertex(reflectionCoefficient * contractionCoefficient, simplexCentroid, worstVertex);
        const FloatType outsideContractedValue = function(outsideContractedVertex);
        if(outsideContractedValue <= reflectedValue) {
          replaceWorst(values, outsideContractedVertex, outsideContractedValue, vertices);
        } else {
          shrink(vertices, values, function);
        }
      } else {
        /* Inside contraction */
        const VectorType insideContractedVertex = generateVertex(-contractionCoefficient, simplexCentroid, worstVertex);
        const FloatType insideContractedValue = function(insideContractedVertex);
        if(insideContractedValue < worstVertexValue) {
          replaceWorst(values, insideContractedVertex, insideContractedValue, vertices);
        } else {
          shrink(vertices, values, function);
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

} // namespace temple
} // namespace Scine

#endif
