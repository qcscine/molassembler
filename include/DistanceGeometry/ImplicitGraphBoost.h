#ifndef INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_TRAITS_H
#define INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_TRAITS_H

#include "DistanceGeometry/ImplicitGraphHelpers.h"
#include "boost/graph/graph_traits.hpp"

/*! @file
 *
 * Includes graph trait, property map and property trait type definitions
 * necessary for ImplicitGraph interoperability with boost::graph
 * algorithms.
 */

/* Boost graph inclusion additions */
namespace boost {

/*!
 * Although the BGL documentation says vertex_and_edge_list_graph_tag is a
 * permissible tag, it does not exist in graph/graph_traits.hpp. So we have to
 * invent one that is convertible to vertex_list_graph_tag and edge_list_graph
 * tag here:
 */
struct vertex_and_edge_list_plus_incidence_graph_tag 
  : virtual vertex_list_graph_tag, 
    virtual edge_list_graph_tag,
    virtual incidence_graph_tag {};

template<>
struct graph_traits<const MoleculeManip::DistanceGeometry::ImplicitGraph> {
  // Shortcut typedef
  using T = MoleculeManip::DistanceGeometry::ImplicitGraph;

  using vertex_descriptor = T::VertexDescriptor;
  using edge_descriptor = T::EdgeDescriptor;

  using directed_category = directed_tag;
  using edge_parallel_category = disallow_parallel_edge_tag;
  using traversal_category = vertex_and_edge_list_plus_incidence_graph_tag;

  using vertices_size_type = T::VertexDescriptor;
  using edges_size_type = T::VertexDescriptor;
  using degree_size_type = T::VertexDescriptor;

  static vertex_descriptor null_vertex() {
    return std::numeric_limits<vertex_descriptor>::max();
  }

  using vertex_iterator = T::vertex_iterator;
  using edge_iterator = T::edge_iterator;
  using out_edge_iterator = T::edge_iterator;
};

template<>
struct graph_traits<MoleculeManip::DistanceGeometry::ImplicitGraph> {
  // Shortcut typedef
  using T = MoleculeManip::DistanceGeometry::ImplicitGraph;

  using vertex_descriptor = T::VertexDescriptor;
  using edge_descriptor = T::EdgeDescriptor;

  using directed_category = directed_tag;
  using edge_parallel_category = disallow_parallel_edge_tag;
  using traversal_category = vertex_and_edge_list_plus_incidence_graph_tag;

  using vertices_size_type = T::VertexDescriptor;
  using edges_size_type = T::VertexDescriptor;
  using degree_size_type = T::VertexDescriptor;

  static vertex_descriptor null_vertex() {
    return std::numeric_limits<vertex_descriptor>::max();
  }

  using vertex_iterator = T::vertex_iterator;
  using edge_iterator = T::edge_iterator;
  using out_edge_iterator = T::edge_iterator;
};

/* edge_weight_t property_map */

template<>
struct property_map<MoleculeManip::DistanceGeometry::ImplicitGraph, edge_weight_t> {
  using T = MoleculeManip::DistanceGeometry::ImplicitGraph;

  using type = T::EdgeWeightMap;
  using const_type = const T::EdgeWeightMap;
};

template<>
struct property_traits<MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeWeightMap> {
  using T = MoleculeManip::DistanceGeometry::ImplicitGraph;

  using value_type = double;
  using reference = double;
  using key_type = T::EdgeDescriptor;
  using category = read_write_property_map_tag;
};

template<>
struct property_traits<const MoleculeManip::DistanceGeometry::ImplicitGraph::EdgeWeightMap> {
  using T = MoleculeManip::DistanceGeometry::ImplicitGraph;

  using value_type = double;
  using reference = double;
  using key_type = T::EdgeDescriptor;
  using category = readable_property_map_tag;
};

/* vertex_index_t property map */

template<>
struct property_map<MoleculeManip::DistanceGeometry::ImplicitGraph, vertex_index_t> {
  using T = MoleculeManip::DistanceGeometry::ImplicitGraph;

  using type = T::VertexIndexMap;
  using const_type = const T::VertexIndexMap;
};

template<>
struct property_traits<MoleculeManip::DistanceGeometry::ImplicitGraph::VertexIndexMap> {
  using T = MoleculeManip::DistanceGeometry::ImplicitGraph;

  using value_type = graph_traits<T>::vertex_descriptor;
  using reference = graph_traits<T>::vertex_descriptor;
  using key_type = graph_traits<T>::vertex_descriptor;
  using category = readable_property_map_tag;
};

template<>
struct property_traits<const MoleculeManip::DistanceGeometry::ImplicitGraph::VertexIndexMap> {
  using T = MoleculeManip::DistanceGeometry::ImplicitGraph;

  using value_type = graph_traits<T>::vertex_descriptor;
  using reference = graph_traits<T>::vertex_descriptor;
  using key_type = graph_traits<T>::vertex_descriptor;
  using category = readable_property_map_tag;
};

} // namespace boost

#endif
