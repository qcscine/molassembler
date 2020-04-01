/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Algorithms to help in the interpretation of the molecular graph
 */

#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H

#include <vector>
#include "Types.h"

namespace Scine {
namespace molassembler {

// Forward-declarations
class Graph;

/*! @brief Calculates the graph distance from a single atom index to all others
 *
 * Performs a BFS through the entire graph starting at the supplied index and
 * records the distance of vertices it encounters.
 *
 * @complexity{@math{O(N)}}
 *
 * @throws std::out_of_range If i >= N()
 *
 * @returns A vector containing the distances of all vertices to the supplied
 *   index
 */
MASM_EXPORT std::vector<unsigned> distance(AtomIndex i, const Graph& graph);

} // namespace molassembler
} // namespace Scine

#endif
