/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Algorithms to help in the interpretation of the molecular graph
 */

#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H

#include "Utils/Geometry/ElementTypes.h"
#include "Types.h"
#include "boost/functional/hash.hpp"

#include <unordered_map>
#include <vector>
#include <memory>
#include <limits>

namespace Scine {
namespace Molassembler {

// Forward-declarations
class Graph;

/*! @brief Calculates the graph distance from a single atom index to all others
 *
 * Performs a BFS through the entire graph starting at the supplied index and
 * records the distance of vertices it encounters.
 *
 * @complexity{@math{O(N)}}
 *
 * @throws std::out_of_range If source >= N()
 *
 * @returns A vector containing the distances of all vertices to the supplied
 *   index
 */
MASM_EXPORT std::vector<unsigned> distance(AtomIndex source, const Graph& graph);

struct MASM_EXPORT PredecessorMap {
  std::vector<AtomIndex> predecessors;

  /**
   * @brief Generate path from source to target vertex
   *
   * @param target Target vertex of shortest path
   *
   * @return Path starting at source and ending at target, including both
   *   source and target vertices
   */
  std::vector<AtomIndex> path(AtomIndex target) const;
};

/**
 * @brief Generates shortest paths to each vertex in a graph
 *
 * @param source Vertex to start at
 * @param graph Graph containing the vertex
 *
 * @throws std::out_of_range If source >= N()
 *
 * @return A flat predecessor map
 */
MASM_EXPORT PredecessorMap shortestPaths(AtomIndex source, const Graph& graph);

/**
 * @brief Base class for edit cost functions
 */
struct MASM_EXPORT EditCost {
  virtual ~EditCost() = default;

  // Alteration means insertion or deletion
  virtual unsigned vertexAlteration() const = 0;
  virtual unsigned edgeAlteration() const = 0;
  virtual unsigned elementSubstitution(Utils::ElementType from, Utils::ElementType to) const = 0;
  virtual unsigned bondSubstitution(BondType from, BondType to) const = 0;
};

/**
 * @brief Fuzzy-matching cost function that prefers element substitutions to
 *   alterations
 */
struct MASM_EXPORT FuzzyCost : public EditCost {
  unsigned vertexAlteration() const override { return 1; }
  unsigned edgeAlteration() const override { return 1; }
  unsigned elementSubstitution(Utils::ElementType from, Utils::ElementType to) const override {
    return (from != to) ? 1 : 0;
  }
  unsigned bondSubstitution(BondType from, BondType to) const override {
    return (from != to) ? 1 : 0;
  }
};

/**
 * @brief Element type conserving cost function, really allowing only edge
 *   alterations
 */
struct MASM_EXPORT ElementsConservedCost : public FuzzyCost {
  unsigned vertexAlteration() const override { return 100; }
  unsigned elementSubstitution(Utils::ElementType from, Utils::ElementType to) const override {
    return (from != to) ? 100 : 0;
  }
};

/**
 * @brief Result type for single-pair graph edit distance calculation
 */
struct MASM_EXPORT Edits {
  using IndexMap = std::vector<AtomIndex>;
  //! Sentinel value for non-existent vertex
  static constexpr AtomIndex epsilon = std::numeric_limits<AtomIndex>::max();

  //! Total distance between graphs according to cost function
  unsigned distance;
  IndexMap indexMap;
};

/**
 * @brief Exact graph edit distance calculation
 *
 * Graph edit distance is symmetric, so order of arguments is irrelevant.
 *
 * @param a First graph to calculate edit distance for
 * @param b Second graph to calculate edit distance for
 * @param cost Cost function implementation for distance calculation
 *
 * @return Graph edit distance metric
 */
Edits minimalEdits(const Graph& a, const Graph& b, std::unique_ptr<EditCost> cost = std::make_unique<FuzzyCost>());

struct MultiEdits {
  using ComponentIndexPair = std::pair<unsigned, AtomIndex>;
  unsigned distance;

  std::unordered_map<ComponentIndexPair, ComponentIndexPair, boost::hash<ComponentIndexPair>> indexMap;
};

using GraphList = std::vector<std::reference_wrapper<Graph>>;
MultiEdits reactionEdits(const GraphList& lhs, const GraphList& rhs);

} // namespace Molassembler
} // namespace Scine

#endif
