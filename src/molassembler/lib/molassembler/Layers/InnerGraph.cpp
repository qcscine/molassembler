#include "molassembler/Layers/Inner.h"

namespace molassembler {

BondType InnerGraph::bondType(const InnerGraph::Edge& edge) const {
  return _graph[edge].bondType;
}

Delib::ElementType InnerGraph::elementType(const Vertex a) const {
  return _graph[a].elementType;
}

InnerGraph::Edge InnerGraph::edge(const Vertex a, const Vertex b) const {
  auto edge = boost::edge(a, b, _graph);
  assert(edge.second);
  return edge.first;
}

InnerGraph::Vertex InnerGraph::source(const InnerGraph::Edge& edge) const {
  return boost::source(edge, _graph);
}

InnerGraph::Vertex InnerGraph::target(const InnerGraph::Edge& edge) const {
  return boost::target(edge, _graph);
}

InnerGraph::VertexRange InnerGraph::vertices() const {
  return boost::vertices(_graph);
}

InnerGraph::EdgeRange InnerGraph::edges() const {
  return boost::edges(_graph);
}

InnerGraph::AdjacentVertexRange InnerGraph::adjacents(const Vertex a) const {
  return boost::adjacent_vertices(a, _graph);
}

InnerGraph::IncidentEdgeRange InnerGraph::edges(const Vertex a) const {
  return boost::out_edges(a, _graph);
}

InnerGraph::BGLType& InnerGraph::bgl() {
  return _graph;
}

const InnerGraph::BGLType& InnerGraph::bgl() const {
  return _graph;
}

} // namespace molassembler
