#ifndef INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H
#define INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H

#include "molassembler/Layers/Outer.h"
#include "molassembler/Layers/Inner.h"

namespace molassembler {

inline InnerGraph::Edge toInner(const OuterGraph::BondIndex& bondIndex, const InnerGraph& graph) {
  return graph.edge(bondIndex.first, bondIndex.second);
}

inline OuterGraph::BondIndex toOuter(const InnerGraph::Edge& edge, const InnerGraph& graph) {
  return { graph.source(edge), graph.target(edge) };
}

} // namespace molassembler

#endif
