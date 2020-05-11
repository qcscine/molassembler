/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Traits definitions for use of ImplicitBoundsGraph with Boost Graph Library
 *
 * Includes graph trait, property map and property trait type definitions
 * necessary for ImplicitBoundsGraph interoperability with boost::graph
 * algorithms.
 */

#ifndef INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_TRAITS_H
#define INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_GRAPH_BOOST_TRAITS_H

// DO NOT CHANGE THIS INCLUDE ORDER
#include "molassembler/DistanceGeometry/ImplicitBoundsGraphHelpers.h"
#include "boost/graph/graph_traits.hpp"

/**
 * @brief Additions to Boost namespace for BGL compatibility of custom graph
 *   classes
 */
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
struct graph_traits<const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph> {
  // Shortcut typedef
  using T = Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph;

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
struct graph_traits<Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph> {
  // Shortcut typedef
  using T = Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph;

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
struct property_map<Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph, edge_weight_t> {
  using T = Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph;

  using type = T::EdgeWeightMap;
  using const_type = const T::EdgeWeightMap;
};

template<>
struct property_traits<Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeWeightMap> {
  using T = Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph;

  using value_type = double;
  using reference = double;
  using key_type = T::EdgeDescriptor;
  using category = read_write_property_map_tag;
};

template<>
struct property_traits<const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::EdgeWeightMap> {
  using T = Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph;

  using value_type = double;
  using reference = double;
  using key_type = T::EdgeDescriptor;
  using category = readable_property_map_tag;
};

/* vertex_index_t property map */

template<>
struct property_map<Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph, vertex_index_t> {
  using T = Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph;

  using type = T::VertexIndexMap;
  using const_type = const T::VertexIndexMap;
};

template<>
struct property_traits<Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexIndexMap> {
  using T = Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph;

  using value_type = graph_traits<T>::vertex_descriptor;
  using reference = graph_traits<T>::vertex_descriptor;
  using key_type = graph_traits<T>::vertex_descriptor;
  using category = readable_property_map_tag;
};

template<>
struct property_traits<const Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph::VertexIndexMap> {
  using T = Scine::Molassembler::DistanceGeometry::ImplicitBoundsGraph;

  using value_type = graph_traits<T>::vertex_descriptor;
  using reference = graph_traits<T>::vertex_descriptor;
  using key_type = graph_traits<T>::vertex_descriptor;
  using category = readable_property_map_tag;
};

} // namespace boost

#endif
