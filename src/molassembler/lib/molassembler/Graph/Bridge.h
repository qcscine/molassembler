/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief OuterGraph and InnerGraph vertex- & edge descriptor conversions
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H
#define INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H

#include "molassembler/OuterGraph.h"
#include "molassembler/Graph/InnerGraph.h"

namespace Scine {

namespace molassembler {

inline InnerGraph::Edge toInner(const BondIndex& bondIndex, const InnerGraph& graph) {
  return graph.edge(bondIndex.first, bondIndex.second);
}

inline BondIndex toOuter(const InnerGraph::Edge& edge, const InnerGraph& graph) {
  return { graph.source(edge), graph.target(edge) };
}

} // namespace molassembler

} // namespace Scine

#endif
