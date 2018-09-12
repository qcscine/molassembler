#ifndef INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H
#define INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H

#include "molassembler/OuterGraph.h"
#include "molassembler/Graph/InnerGraph.h"

/*!@file
 *
 * @brief OuterGraph and InnerGraph vertex- & edge descriptor conversions
 */

namespace molassembler {

inline InnerGraph::Edge toInner(const BondIndex& bondIndex, const InnerGraph& graph) {
  return graph.edge(bondIndex.first, bondIndex.second);
}

inline BondIndex toOuter(const InnerGraph::Edge& edge, const InnerGraph& graph) {
  return { graph.source(edge), graph.target(edge) };
}

} // namespace molassembler

#endif
