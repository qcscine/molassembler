/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief GOR1 single-source shortest-paths algorithm with boost graph
 *
 * Implements a simplified GOR1 single-source shortest paths algorithm in the
 * boost namespace in Boost Graph Library style (including a Visitor for
 * algorithm observation and side effects).
 *
 * Reference:
 * - Cherkassky, B. V., Goldberg, A. V., & Radzik, T. (1996). Shortest paths
 *   algorithms: Theory and experimental evaluation. Mathematical Programming,
 *   73(2), 129–174. https://doi.org/10.1007/BF02592101
 */

#include <stack>
#include <limits>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/graph_concepts.hpp"

namespace boost {

namespace detail {

/*!
 * @brief A Dummy Visitor object demonstrating the interface needed to enable
 *   GOR1 algorithm observation and side effects.
 */
struct DummyGor1Visitor {
/* Stack operations */
  //! Dummy stack A push event caller
  template<typename VertexDescriptor, typename IncidenceGraph>
  void a_push(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  //! Dummy stack A pop event caller
  template<typename VertexDescriptor, typename IncidenceGraph>
  void a_pop(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  //! Dummy stack A push event caller
  template<typename VertexDescriptor, typename IncidenceGraph>
  void b_push(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  //! Dummy stack A pop event caller
  template<typename VertexDescriptor, typename IncidenceGraph>
  void b_pop(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

/* Edge operations */
  //! Dummy edge examine event caller
  template<typename EdgeDescriptor, typename IncidenceGraph>
  void examine_edge(const EdgeDescriptor& /* e */, const IncidenceGraph& /* g */) {}

  //! Dummy edge relax event caller
  template<typename EdgeDescriptor, typename IncidenceGraph>
  void relax_edge(const EdgeDescriptor& /* e */, const IncidenceGraph& /* g */) {}

/* Vertex colorings */
  //! Dummy vertex marked white event caller
  template<typename VertexDescriptor, typename IncidenceGraph>
  void mark_white(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  //! Dummy vertex marked gray event caller
  template<typename VertexDescriptor, typename IncidenceGraph>
  void mark_gray(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  //! Dummy vertex marked blac event caller
  template<typename VertexDescriptor, typename IncidenceGraph>
  void mark_black(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}
};

/*!
 * @brief GOR1 helper function that performs the scanning of a vertex' edges
 *
 * @tparam VertexDescriptor Type of the Graph's vertex descriptor
 * @tparam IncidenceGraph Type modeling Boost's IncidenceGraph concept
 * @tparam DistanceMap Type mapping vertex descriptors to distance
 * @tparam PredecessorMap Type mapping vertex descriptors to their predecessor
 *   vertex
 * @tparam ColorMap Type mapping vertex descriptors to a color
 * @tparam Visitor Type that performs event visitation operations
 *
 * @param vertex The vertex to scan
 * @param graph The graph that vertex is contained in
 * @param predecessor_map The map of vertices to their predecessors
 * @param color_map The map of vertices to their color
 * @param distance_map The map of vertices to their distance
 * @param B The B stack as described in the paper
 * @param visitor A visitor that matches the DummyGor1Visitor interface
 */
template<
  typename VertexDescriptor,
  class IncidenceGraph,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap,
  class Visitor
>
void gor1_simplified_scan(
  const VertexDescriptor& vertex,
  const IncidenceGraph& graph,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map,
  std::stack<VertexDescriptor>& B,
  Visitor& visitor
) {
  using ColorValue = typename property_traits<ColorMap>::value_type;
  using Color = color_traits<ColorValue>;

  auto vertexDistance = get(distance_map, vertex);
  auto out_iter_pair = out_edges(vertex, graph);

  while(out_iter_pair.first != out_iter_pair.second) {
    auto edge = *out_iter_pair.first;
    visitor.examine_edge(edge, graph);

    auto targetVertex = target(edge, graph);
    auto edgeWeight = get(boost::edge_weight, graph, edge);

    if(vertexDistance + edgeWeight < get(distance_map, targetVertex)) {
      // Set distance for target vertex
      put(distance_map, targetVertex, vertexDistance + edgeWeight);

      // Mark target as vertex's descendant
      put(predecessor_map, targetVertex, vertex);

      visitor.relax_edge(edge, graph);

      // If not marked scanned, add it to B and mark it
      if(get(color_map, targetVertex) != Color::black()) {
        B.push(targetVertex);
        visitor.b_push(targetVertex, graph);
        put(color_map, targetVertex, Color::gray());
        visitor.mark_gray(targetVertex, graph);
      }
    }

    ++out_iter_pair.first;
  }
}

} // namespace detail

/*!
 * @brief Simplified GOR1 single source shortest paths algorithm
 *
 * Implements the algorithm described in
 * - Cherkassky, B. V., Goldberg, A. V., & Radzik, T. (1996). Shortest paths
 *   algorithms: Theory and experimental evaluation. Mathematical Programming,
 *   73(2), 129–174. https://doi.org/10.1007/BF02592101
 *
 * Reference C implementation in
 * - https://github.com/skvadrik/cherkassky_goldberg_radzik
 *
 * @complexity{@math{\Theta(V E)}}
 *
 * @tparam IncidenceGraph Type modeling Boost's IncidenceGraph concept
 * @tparam DistanceMap Type mapping vertex descriptors to distance
 * @tparam PredecessorMap Type mapping vertex descriptors to their predecessor
 *   vertex
 * @tparam ColorMap Type mapping vertex descriptors to a color
 * @tparam Visitor Type that performs event visitation operations
 * @tparam VertexDescriptor Type of the Graph's vertex descriptor
 *
 * @param graph The graph that vertex is contained in
 * @param root_vertex The source vertex from which shortest distances are to be
 *   calculated
 * @param predecessor_map The map of vertices to their predecessors
 * @param color_map The map of vertices to their color
 * @param distance_map The map of vertices to their distance
 * @param visitor A visitor that matches the DummyGor1Visitor interface
 *
 * @note This function, in boost style, expects you to set up some of the
 * data structures used internally in the algorithm, whether you care about
 * them or not for the off chance you want all of them and structured bindings
 * weren't a thing when boost came about.
 *
 * @code{.cpp}
 * // Setup prior to call
 * const unsigned N = boost::num_vertices(g);
 * std::vector<double> distances(N);
 * std::vector<VertexType> predecessors(N);
 *
 * auto predecessor_map = boost::make_iterator_property_map(
 *   predecessors.begin(),
 *   get(boost::vertex_index, g)
 * );

 * auto distance_map = boost::make_iterator_property_map(
 *   distances.begin(),
 *   get(boost::vertex_index, g)
 * );

 * using ColorMapType = boost::two_bit_color_map<>;
 * ColorMapType color_map {N};
 * @endcode
 */
template<
  class IncidenceGraph,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap,
  class Visitor,
  typename VertexDescriptor
>
bool gor1_simplified_shortest_paths(
  const IncidenceGraph& graph,
  const VertexDescriptor& root_vertex,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map,
  Visitor& visitor
) {
  BOOST_CONCEPT_ASSERT(( IncidenceGraphConcept<IncidenceGraph> ));

  using Vertex = typename graph_traits<IncidenceGraph>::vertex_descriptor;
  BOOST_CONCEPT_ASSERT(( ReadWritePropertyMapConcept<ColorMap, Vertex> ));

  using Edge = typename graph_traits<IncidenceGraph>::edge_descriptor;
  using EdgeMap = typename property_map<IncidenceGraph, edge_weight_t>::type;
  BOOST_CONCEPT_ASSERT(( ReadablePropertyMapConcept<EdgeMap, Edge> ));

  using ColorValue = typename property_traits<ColorMap>::value_type;
  using Color = color_traits<ColorValue>;

  auto verticesIterPair = vertices(graph);
  while(verticesIterPair.first != verticesIterPair.second) {
    put(distance_map, *verticesIterPair.first, std::numeric_limits<double>::max());
    ++verticesIterPair.first;
  }

  put(distance_map, root_vertex, 0.0);
  put(color_map, root_vertex, Color::gray());
  put(predecessor_map, root_vertex, root_vertex);

  std::stack<VertexDescriptor> A, B;
  B.push(root_vertex);
  visitor.b_push(root_vertex, graph);

  while(!B.empty()) {
    // Compute A from B, emptying B in the process
    while(!B.empty()) {
      VertexDescriptor v = B.top();
      B.pop();

      visitor.b_pop(v, graph);

      ColorValue v_color = get(color_map, v);

      if(v_color == Color::black()) {
        A.push(v);
        visitor.a_push(v, graph);
      } else if(v_color == Color::gray()) {
        B.push(v);
        visitor.b_push(v, graph);
        put(color_map, v, Color::black());
        visitor.mark_black(v, graph);

        detail::gor1_simplified_scan(
          v,
          graph,
          predecessor_map,
          color_map,
          distance_map,
          B,
          visitor
        );
      }
    }

    // Scan all elements in A, populating B with nodes added to the tree
    while(!A.empty()) {
      VertexDescriptor v = A.top();
      A.pop();

      visitor.a_pop(v, graph);

      // Scan
      detail::gor1_simplified_scan(
        v,
        graph,
        predecessor_map,
        color_map,
        distance_map,
        B,
        visitor
      );

      // Mark white
      put(color_map, v, Color::white());
      visitor.mark_white(v, graph);
    }
  }

  return true;
}

//! @overload
template<
  class IncidenceGraph,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap,
  typename VertexDescriptor
>
bool gor1_simplified_shortest_paths(
  const IncidenceGraph& graph,
  const VertexDescriptor& root_vertex,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map
) {
  /* This function is needed since binding an rvalue to a lvalue reference is
   * not permitted in:
   *
   * template<..., Visitor = detail::DummyGor1Visitor>
   * bool shortest_paths(
   *   ...,
   *   Visitor& visitor = Visitor {}
   * )
   */
  detail::DummyGor1Visitor visitor;

  return gor1_simplified_shortest_paths(
    graph,
    root_vertex,
    predecessor_map,
    color_map,
    distance_map,
    visitor
  );
}


} // namespace boost
