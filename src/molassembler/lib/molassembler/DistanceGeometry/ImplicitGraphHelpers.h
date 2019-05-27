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

// DO NOT CHANGE THIS INCLUDE ORDER
#include "molassembler/DistanceGeometry/ImplicitGraph.h"
#include "boost/graph/properties.hpp"

namespace boost {

/* VertexListGraph concept */
inline std::pair<
  Scine::molassembler::DistanceGeometry::ImplicitGraph::vertex_iterator,
  Scine::molassembler::DistanceGeometry::ImplicitGraph::vertex_iterator
> vertices(
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return {
    g.vbegin(),
    g.vend()
  };
}

inline Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor num_vertices(
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.num_vertices();
}

/* EdgeListGraph concept */
inline std::pair<
  Scine::molassembler::DistanceGeometry::ImplicitGraph::edge_iterator,
  Scine::molassembler::DistanceGeometry::ImplicitGraph::edge_iterator
> edges(const Scine::molassembler::DistanceGeometry::ImplicitGraph& g) {
  return {
    g.ebegin(),
    g.eend()
  };
}

inline Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor num_edges(
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.num_edges();
}

inline Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor source(
  const Scine::molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor& e,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.source(e);
}

inline Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor target(
  const Scine::molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor& e,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.target(e);
}

/* AdjacencyMatrix concept */
inline std::pair<
  Scine::molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor,
  bool
> edge(
  const Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor& u,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor& v,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.edge(u, v);
}

/* PropertyGraph concept */
inline const Scine::molassembler::DistanceGeometry::ImplicitGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline Scine::molassembler::DistanceGeometry::ImplicitGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline double get(
  const boost::edge_weight_t& p,
  Scine::molassembler::DistanceGeometry::ImplicitGraph& g,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline double get(
  const boost::edge_weight_t& p,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline void put(
  const boost::edge_weight_t&,
  Scine::molassembler::DistanceGeometry::ImplicitGraph&,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph::EdgeDescriptor&,
  double
) {
  /* do nothing */
}

inline Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexIndexMap get(
  const boost::vertex_index_t&,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph&
) {
  return Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexIndexMap {};
}

/* IncidenceGraph concept */
inline std::pair<
  Scine::molassembler::DistanceGeometry::ImplicitGraph::edge_iterator,
  Scine::molassembler::DistanceGeometry::ImplicitGraph::edge_iterator
> out_edges(
  const Scine::molassembler::DistanceGeometry::ImplicitGraph::VertexDescriptor& u,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return {
    g.obegin(u),
    g.oend(u)
  };
}

inline unsigned long out_degree(
  const unsigned long& u,
  const Scine::molassembler::DistanceGeometry::ImplicitGraph& g
) {
  return g.out_degree(u);
}

} // namespace boost

#endif
