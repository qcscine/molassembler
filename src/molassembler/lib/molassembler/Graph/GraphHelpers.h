#ifndef INCLUDE_MOLASSEMBLER_GRAPH_HELPERS_H
#define INCLUDE_MOLASSEMBLER_GRAPH_HELPERS_H

#include "molassembler/Detail/RangeForTemporary.h"
#include "molassembler/Types.h"

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
inline AtomIndex addAtom(
  const Delib::ElementType elementType,
  OuterGraph& graph
) {
  auto vertex = boost::add_vertex(graph);
  graph[vertex].elementType = elementType;
  return vertex;
}

inline BondIndex addBond(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType,
  OuterGraph& graph
) {
  auto edgeAddPair = boost::add_edge(a, b, graph);

  if(!edgeAddPair.second) {
    throw std::logic_error("Cannot add a bond where one already is present!");
  }

  graph[edgeAddPair.first].bondType = bondType;

  return edgeAddPair.first;
}

inline RangeForTemporary<OuterGraph::vertex_iterator> vertices(const OuterGraph& graph) {
  return RangeForTemporary<OuterGraph::vertex_iterator>(
    boost::vertices(graph)
  );
}

inline RangeForTemporary<OuterGraph::edge_iterator> edges(const OuterGraph& graph) {
  return RangeForTemporary<OuterGraph::edge_iterator>(
    boost::edges(graph)
  );
}

inline RangeForTemporary<OuterGraph::out_edge_iterator> edges(
  const AtomIndex a,
  const OuterGraph& graph
) {
  return RangeForTemporary<OuterGraph::out_edge_iterator>(
    boost::out_edges(a, graph)
  );
}

inline RangeForTemporary<OuterGraph::adjacency_iterator> adjacents(
  const AtomIndex a,
  const OuterGraph& graph
) {
  return RangeForTemporary<OuterGraph::adjacency_iterator>(
    boost::adjacent_vertices(a, graph)
  );
}

inline AtomIndex source(
  const BondIndex edge,
  const OuterGraph& graph
) {
  return boost::source(edge, graph);
}

inline AtomIndex target(
  const BondIndex edge,
  const OuterGraph& graph
) {
  return boost::target(edge, graph);
}

inline Delib::ElementType& elementType(const AtomIndex a, OuterGraph& graph) {
  return graph[a].elementType;
}

inline Delib::ElementType elementType(const AtomIndex a, const OuterGraph& graph) {
  return graph[a].elementType;
}

inline BondIndex edge(const AtomIndex a, const AtomIndex b, const OuterGraph& graph) {
  auto edgePair = boost::edge(a, b, graph);
  assert(edgePair.second);
  return edgePair.first;
}

inline boost::optional<BondIndex> edgeOption(
  const AtomIndex a,
  const AtomIndex b,
  const OuterGraph& graph
) {
  auto edgePair = boost::edge(a, b, graph);

  if(edgePair.second) {
    return edgePair.first;
  }

  // fallback
  return boost::none;
}

inline BondType& bondType(const BondIndex edge, OuterGraph& graph) {
  return graph[edge].bondType;
}

inline BondType bondType(const BondIndex edge, const OuterGraph& graph) {
  return graph[edge].bondType;
}

inline unsigned numAdjacencies(const AtomIndex a, const OuterGraph& graph) {
  return boost::out_degree(a, graph);
}

inline AtomIndex numVertices(const OuterGraph& graph) {
  return boost::num_vertices(graph);
}

} // namespace graph

} // namespace molassembler

#endif
