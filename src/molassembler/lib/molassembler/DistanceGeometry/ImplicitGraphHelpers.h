/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief boost-namespace functions for use of ImplicitGraph with BGL
 *
 * Includes function definitions necessary for ImplicitGraph
 * interoperability with boost::graph algorithms.
 */

#ifndef INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_CONCEPTS_H
#define INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_CONCEPTS_H

#include "molassembler/DistanceGeometry/ImplicitGraph.h"
#include "boost/graph/properties.hpp"

namespace boost {

/* VertexListGraph concept */
inline std::pair<
  molassembler::DistanceGeometry::ImplicitGraph::vertex_iterator,
  molassembler::DistanceGeometry::ImplicitGraph::vertex_iterator
> vertices(
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return {
    g.vbegin(),
    g.vend()
  };
}

inline molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor num_vertices(
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.num_vertices();
}

/* EdgeListGraph concept */
inline std::pair<
  molassembler::DistanceGeometry::ImplicitGraph::edge_iterator,
  molassembler::DistanceGeometry::ImplicitGraph::edge_iterator
> edges(const molassembler::DistanceGeometry::ImplicitGraph& g) {
  return {
    g.ebegin(),
    g.eend()
  };
}

inline molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor num_edges(
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.num_edges();
}

inline molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor source(
  const molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor& e,
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.source(e);
}

inline molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor target(
  const molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor& e,
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.target(e);
}

/* AdjacencyMatrix concept */
inline std::pair<
  molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor,
  bool
> edge(
  const molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor& u,
  const molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor& v,
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.edge(u, v);
}

/* PropertyGraph concept */
inline const molassembler::DistanceGeometry::ImplicitGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline molassembler::DistanceGeometry::ImplicitGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline double get(
  const boost::edge_weight_t& p,
  molassembler::DistanceGeometry::ImplicitGraph& g,
  const molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline double get(
  const boost::edge_weight_t& p,
  const molassembler::DistanceGeometry::ImplicitGraph& g,
  const molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline void put(
  const boost::edge_weight_t&,
  molassembler::DistanceGeometry::ImplicitGraph&,
  const molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor&,
  double
) {
  /* do nothing */
}

inline molassembler::DistanceGeometry::ImplicitGraph::VertexIndexMap get(
  const boost::vertex_index_t&,
  const molassembler::DistanceGeometry::ImplicitGraph&
) {
  return molassembler::DistanceGeometry::ImplicitGraph::VertexIndexMap {};
}

/* IncidenceGraph concept */
inline std::pair<
  molassembler::DistanceGeometry::ImplicitGraph::edge_iterator,
  molassembler::DistanceGeometry::ImplicitGraph::edge_iterator
> out_edges(
  const molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor& u,
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return {
    g.obegin(u),
    g.oend(u)
  };
}

inline unsigned long out_degree(
  const unsigned long& u,
  const molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.out_degree(u);
}

} // namespace boost

#endif
