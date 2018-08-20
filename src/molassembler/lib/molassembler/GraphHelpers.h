#ifndef INCLUDE_MOLASSEMBLER_GRAPH_HELPERS_H
#define INCLUDE_MOLASSEMBLER_GRAPH_HELPERS_H

#include "molassembler/detail/SharedTypes.h"
#include "molassembler/detail/RangeForTemporary.h"

/*!@file
 *
 * Contains common graph manipulation and information stubs
 */

namespace molassembler {

namespace graph {

/*!
 * Adds an vertex to the graph and sets it's element type property, returning
 * the new index.
 */
inline AtomIndexType addAtom(
  const Delib::ElementType elementType,
  GraphType& graph
) {
  auto vertex = boost::add_vertex(graph);
  graph[vertex].elementType = elementType;
  return vertex;
}

inline GraphType::edge_descriptor addBond(
  const AtomIndexType a,
  const AtomIndexType b,
  const BondType bondType,
  GraphType& graph
) {
  auto edgeAddPair = boost::add_edge(a, b, graph);

  if(!edgeAddPair.second) {
    throw std::logic_error("Cannot add a bond where one already is present!");
  }

  graph[edgeAddPair.first].bondType = bondType;

  return edgeAddPair.first;
}

inline RangeForTemporary<GraphType::vertex_iterator> vertices(const GraphType& graph) {
  return RangeForTemporary<GraphType::vertex_iterator>(
    boost::vertices(graph)
  );
}

inline RangeForTemporary<GraphType::edge_iterator> edges(const GraphType& graph) {
  return RangeForTemporary<GraphType::edge_iterator>(
    boost::edges(graph)
  );
}

inline RangeForTemporary<GraphType::out_edge_iterator> edges(
  const AtomIndexType a,
  const GraphType& graph
) {
  return RangeForTemporary<GraphType::out_edge_iterator>(
    boost::out_edges(a, graph)
  );
}

inline RangeForTemporary<GraphType::adjacency_iterator> adjacents(
  const AtomIndexType a,
  const GraphType& graph
) {
  return RangeForTemporary<GraphType::adjacency_iterator>(
    boost::adjacent_vertices(a, graph)
  );
}

inline AtomIndexType source(
  const GraphType::edge_descriptor edge,
  const GraphType& graph
) {
  return boost::source(edge, graph);
}

inline AtomIndexType target(
  const GraphType::edge_descriptor edge,
  const GraphType& graph
) {
  return boost::target(edge, graph);
}

inline Delib::ElementType& elementType(const AtomIndexType a, GraphType& graph) {
  return graph[a].elementType;
}

inline Delib::ElementType elementType(const AtomIndexType a, const GraphType& graph) {
  return graph[a].elementType;
}

inline GraphType::edge_descriptor edge(const AtomIndexType a, const AtomIndexType b, const GraphType& graph) {
  auto edgePair = boost::edge(a, b, graph);
  assert(edgePair.second);
  return edgePair.first;
}

inline boost::optional<GraphType::edge_descriptor> edgeOption(
  const AtomIndexType a,
  const AtomIndexType b,
  const GraphType& graph
) {
  auto edgePair = boost::edge(a, b, graph);

  if(edgePair.second) {
    return edgePair.first;
  }

  // fallback
  return boost::none;
}

inline BondType& bondType(const GraphType::edge_descriptor edge, GraphType& graph) {
  return graph[edge].bondType;
}

inline BondType bondType(const GraphType::edge_descriptor edge, const GraphType& graph) {
  return graph[edge].bondType;
}

inline unsigned numAdjacencies(const AtomIndexType a, const GraphType& graph) {
  return boost::out_degree(a, graph);
}

inline AtomIndexType numVertices(const GraphType& graph) {
  return boost::num_vertices(graph);
}

} // namespace graph

} // namespace molassembler

#endif
