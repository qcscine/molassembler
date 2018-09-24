#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H

#include <vector>
#include "Types.h"

/*!@file
 *
 * @brief Algorithms to help in the interpretation of the molecular graph
 */

namespace molassembler {

// Forward-declarations
class OuterGraph;

/*! Calculates the graph distance from a single atom index to all others
 *
 * Performs a BFS through the entire graph starting at the supplied index and
 * records the distance of vertices it encounters.
 *
 * @throws std::out_of_range If i >= N()
 *
 * @returns A vector containing the distances of all vertices to the supplied
 *   index
 */
std::vector<unsigned> distance(AtomIndex i, const OuterGraph& graph);

} // namespace molassembler

#endif
