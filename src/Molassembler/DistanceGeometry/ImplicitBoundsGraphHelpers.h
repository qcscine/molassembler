/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief boost-namespace functions for use of ImplicitBoundsGraph with BGL
 *
 * Includes function definitions necessary for ImplicitBoundsGraph
 * interoperability with boost::graph algorithms.
 */

#ifndef INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_CONCEPTS_H
#define INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_CONCEPTS_H

// DO NOT CHANGE THIS INCLUDE ORDER (you'll get defined-after compile errors)
#include "Molassembler/DistanceGeometry/ImplicitBoundsGraph.h"
#include "boost/graph/properties.hpp"

namespace boost {

/* VertexListGraph concept */
inline std::pair<
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::vertex_iterator,
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::vertex_iterator
> vertices(
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return {g.vbegin(), g.vend()};
}

inline Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexDescriptor num_vertices(
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return g.num_vertices();
}

/* EdgeListGraph concept */
inline std::pair<
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::edge_iterator,
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::edge_iterator
> edges(const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g) {
  return {g.ebegin(), g.eend()};
}

inline Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexDescriptor num_edges(
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return g.num_edges();
}

inline Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexDescriptor source(
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeDescriptor& e,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return g.source(e);
}

inline Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexDescriptor target(
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeDescriptor& e,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return g.target(e);
}

/* AdjacencyMatrix concept */
inline std::pair<
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeDescriptor,
  bool
> edge(
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexDescriptor& u,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexDescriptor& v,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return g.edge(u, v);
}

/* PropertyGraph concept */
inline Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeWeightMap get(
  const boost::edge_weight_t& /* tag */,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeWeightMap get(
  const boost::edge_weight_t& /* tag */,
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline double get(
  const boost::edge_weight_t& p,
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline double get(
  const boost::edge_weight_t& p,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline void put(
  const boost::edge_weight_t& /* tag */,
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& /* g */,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeDescriptor& /* e */,
  double /* v */
) {
  /* do nothing */
}

inline Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexIndexMap get(
  const boost::vertex_index_t& /* tag */,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& /* g */
) {
  return Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexIndexMap {};
}

/* IncidenceGraph concept */
inline std::pair<
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::edge_iterator,
  Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::edge_iterator
> out_edges(
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexDescriptor& u,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return {g.obegin(u), g.oend(u)};
}

inline unsigned long out_degree(
  const unsigned long& u,
  const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph& g
) {
  return g.out_degree(u);
}

} // namespace boost

#endif
