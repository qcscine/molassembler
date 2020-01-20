/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Graph and PrivateGraph vertex- & edge descriptor conversions
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H
#define INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H

#include "molassembler/Graph.h"
#include "molassembler/Graph/PrivateGraph.h"

namespace Scine {

namespace molassembler {

//! Transform BondIndex to PrivateGraph::Edge
inline PrivateGraph::Edge toInner(const BondIndex& bondIndex, const PrivateGraph& graph) {
  return graph.edge(bondIndex.first, bondIndex.second);
}

//! Transform PrivateGraph::Edge to BondIndex
inline BondIndex toOuter(const PrivateGraph::Edge& edge, const PrivateGraph& graph) {
  return { graph.source(edge), graph.target(edge) };
}

} // namespace molassembler

} // namespace Scine

#endif
