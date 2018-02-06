#ifndef INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_CONCEPTS_H
#define INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_CONCEPTS_H

/*! @file
 *
 * Includes function definitions necessary for ImplicitGraph
 * interoperability with boost::graph algorithms.
 */

#include "DistanceGeometry/ImplicitGraph.h"
#include "boost/graph/properties.hpp"

namespace boost {

/* VertexListGraph concept */
inline std::pair<
  MoleculeManip::DistanceGeometry::ImplicitGraph::vertex_iterator,
  MoleculeManip::DistanceGeometry::ImplicitGraph::vertex_iterator
> vertices(
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return {
    g.vbegin(),
    g.vend()
  };
}

inline MoleculeManip::DistanceGeometry::ImplicitGraph::VertexDescriptor num_vertices(
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return g.num_vertices();
}

/* EdgeListGraph concept */
inline std::pair<
  MoleculeManip::DistanceGeometry::ImplicitGraph::edge_iterator,
  MoleculeManip::DistanceGeometry::ImplicitGraph::edge_iterator
> edges(const MoleculeManip::DistanceGeometry::ImplicitGraph& g) {
  return {
    g.ebegin(),
    g.eend()
  };
}

inline MoleculeManip::DistanceGeometry::ImplicitGraph::VertexDescriptor num_edges(
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return g.num_edges();
}

inline MoleculeManip::DistanceGeometry::ImplicitGraph::VertexDescriptor source(
  const MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeDescriptor& e,
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return g.source(e);
}

inline MoleculeManip::DistanceGeometry::ImplicitGraph::VertexDescriptor target(
  const MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeDescriptor& e,
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return g.target(e);
}

/* AdjacencyMatrix concept */
inline std::pair<
  MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeDescriptor,
  bool
> edge(
  const MoleculeManip::DistanceGeometry::ImplicitGraph::VertexDescriptor& u,
  const MoleculeManip::DistanceGeometry::ImplicitGraph::VertexDescriptor& v,
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return g.edge(u, v);
}

/* PropertyGraph concept */
inline const MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeWeightMap get(
  const boost::edge_weight_t&,
  MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return g.getEdgeWeightPropertyMap();
}

inline double get(
  const boost::edge_weight_t& p,
  MoleculeManip::DistanceGeometry::ImplicitGraph& g,
  const MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
} 

inline double get(
  const boost::edge_weight_t& p,
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g,
  const MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeDescriptor& x
) {
  return get(p, g)[x];
}

inline void put(
  const boost::edge_weight_t&,
  MoleculeManip::DistanceGeometry::ImplicitGraph&,
  const MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeDescriptor&,
  double
) {
  /* do nothing */
}

inline MoleculeManip::DistanceGeometry::ImplicitGraph::VertexIndexMap get(
  const boost::vertex_index_t&,
  const MoleculeManip::DistanceGeometry::ImplicitGraph&
) {
  return MoleculeManip::DistanceGeometry::ImplicitGraph::VertexIndexMap {};
}

/* IncidenceGraph concept */
inline std::pair<
  MoleculeManip::DistanceGeometry::ImplicitGraph::edge_iterator,
  MoleculeManip::DistanceGeometry::ImplicitGraph::edge_iterator
> out_edges(
  const MoleculeManip::DistanceGeometry::ImplicitGraph::VertexDescriptor& u,
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return {
    g.obegin(u),
    g.oend(u)
  };
}

inline unsigned long out_degree(
  const unsigned long& u,
  const MoleculeManip::DistanceGeometry::ImplicitGraph& g
) {
  return g.out_degree(u);
}

} // namespace boost

#endif
