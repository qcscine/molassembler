/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Algorithms to help in the interpretation of the molecular graph
 */

#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H

#include <vector>
#include "Types.h"

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

struct PredecessorMap {
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

} // namespace Molassembler
} // namespace Scine

#endif
