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

// DO NOT CHANGE THIS INCLUDE ORDER (you'll get defined-after compile errors)
#include "molassembler/DistanceGeometry/ImplicitGraph.h"
#include "boost/graph/properties.hpp"

namespace boost {

/* VertexListGraph concept */
inline std::pair<
  Scine::molassembler::distance_geometry::ImplicitGraph::vertex_iterator,
  Scine::molassembler::distance_geometry::ImplicitGraph::vertex_iterator
> vertices(
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return {
    g.vbegin(),
    g.vend()
  };
}

inline Scine::molassembler::distance_geometry::ImplicitGraph::VertexDescriptor num_vertices(
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return g.num_vertices();
}

/* EdgeListGraph concept */
inline std::pair<
  Scine::molassembler::distance_geometry::ImplicitGraph::edge_iterator,
  Scine::molassembler::distance_geometry::ImplicitGraph::edge_iterator
> edges(const Scine::molassembler::distance_geometry::ImplicitGraph& g) {
  return {
    g.ebegin(),
    g.eend()
  };
}

inline Scine::molassembler::distance_geometry::ImplicitGraph::VertexDescriptor num_edges(
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return g.num_edges();
}

inline Scine::molassembler::distance_geometry::ImplicitGraph::VertexDescriptor source(
  const Scine::molassembler::distance_geometry::ImplicitGraph::EdgeDescriptor& e,
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return g.source(e);
}

inline Scine::molassembler::distance_geometry::ImplicitGraph::VertexDescriptor target(
  const Scine::molassembler::distance_geometry::ImplicitGraph::EdgeDescriptor& e,
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return g.target(e);
}

/* AdjacencyMatrix concept */
inline std::pair<
  Scine::molassembler::distance_geometry::ImplicitGraph::EdgeDescriptor,
  bool
> edge(
  const Scine::molassembler::distance_geometry::ImplicitGraph::VertexDescriptor& u,
  const Scine::molassembler::distance_geometry::ImplicitGraph::VertexDescriptor& v,
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return g.edge(u, v);
}

/* PropertyGraph concept */
inline const Scine::molassembler::distance_geometry::ImplicitGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline Scine::molassembler::distance_geometry::ImplicitGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline double get(
  const boost::edge_weight_t& p,
  Scine::molassembler::distance_geometry::ImplicitGraph& g,
  const Scine::molassembler::distance_geometry::ImplicitGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline double get(
  const boost::edge_weight_t& p,
  const Scine::molassembler::distance_geometry::ImplicitGraph& g,
  const Scine::molassembler::distance_geometry::ImplicitGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline void put(
  const boost::edge_weight_t&,
  Scine::molassembler::distance_geometry::ImplicitGraph&,
  const Scine::molassembler::distance_geometry::ImplicitGraph::EdgeDescriptor&,
  double
) {
  /* do nothing */
}

inline Scine::molassembler::distance_geometry::ImplicitGraph::VertexIndexMap get(
  const boost::vertex_index_t&,
  const Scine::molassembler::distance_geometry::ImplicitGraph&
) {
  return Scine::molassembler::distance_geometry::ImplicitGraph::VertexIndexMap {};
}

/* IncidenceGraph concept */
inline std::pair<
  Scine::molassembler::distance_geometry::ImplicitGraph::edge_iterator,
  Scine::molassembler::distance_geometry::ImplicitGraph::edge_iterator
> out_edges(
  const Scine::molassembler::distance_geometry::ImplicitGraph::VertexDescriptor& u,
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return {
    g.obegin(u),
    g.oend(u)
  };
}

inline unsigned long out_degree(
  const unsigned long& u,
  const Scine::molassembler::distance_geometry::ImplicitGraph& g
) {
  return g.out_degree(u);
}

} // namespace boost

#endif
