/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief boost-namespace functions for use of ImplicitBoundsGraph with BGL
 *
 * Includes function definitions necessary for ImplicitBoundsGraph
 * interoperability with boost::graph algorithms.
 */

#ifndef INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_CONCEPTS_H
#define INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_CONCEPTS_H

// DO NOT CHANGE THIS INCLUDE ORDER (you'll get defined-after compile errors)
#include "molassembler/DistanceGeometry/ImplicitBoundsGraph.h"
#include "boost/graph/properties.hpp"

namespace boost {

/* VertexListGraph concept */
inline std::pair<
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph::vertex_iterator,
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph::vertex_iterator
> vertices(
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return {
    g.vbegin(),
    g.vend()
  };
}

inline Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexDescriptor num_vertices(
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return g.num_vertices();
}

/* EdgeListGraph concept */
inline std::pair<
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph::edge_iterator,
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph::edge_iterator
> edges(const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g) {
  return {
    g.ebegin(),
    g.eend()
  };
}

inline Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexDescriptor num_edges(
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return g.num_edges();
}

inline Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexDescriptor source(
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::EdgeDescriptor& e,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return g.source(e);
}

inline Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexDescriptor target(
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::EdgeDescriptor& e,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return g.target(e);
}

/* AdjacencyMatrix concept */
inline std::pair<
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph::EdgeDescriptor,
  bool
> edge(
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexDescriptor& u,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexDescriptor& v,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return g.edge(u, v);
}

/* PropertyGraph concept */
inline const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline Scine::molassembler::distance_geometry::ImplicitBoundsGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline double get(
  const boost::edge_weight_t& p,
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline double get(
  const boost::edge_weight_t& p,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline void put(
  const boost::edge_weight_t&,
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph&,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::EdgeDescriptor&,
  double
) {
  /* do nothing */
}

inline Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexIndexMap get(
  const boost::vertex_index_t&,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph&
) {
  return Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexIndexMap {};
}

/* IncidenceGraph concept */
inline std::pair<
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph::edge_iterator,
  Scine::molassembler::distance_geometry::ImplicitBoundsGraph::edge_iterator
> out_edges(
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph::VertexDescriptor& u,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return {
    g.obegin(u),
    g.oend(u)
  };
}

inline unsigned long out_degree(
  const unsigned long& u,
  const Scine::molassembler::distance_geometry::ImplicitBoundsGraph& g
) {
  return g.out_degree(u);
}

} // namespace boost

#endif
