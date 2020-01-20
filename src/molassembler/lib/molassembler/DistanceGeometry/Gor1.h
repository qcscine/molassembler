/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief GOR1 specializations for use with ImplicitBoundsGraph and ExplicitBoundsGraph
 *
 * Contains specializations of the Gor1 algorithm for use with ImplicitBoundsGraph
 * and ExplicitBoundsGraph
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_GOR_SPECIALIZATION_H
#define INCLUDE_MOLASSEMBLER_DG_GOR_SPECIALIZATION_H

#include <stack>
#include <limits>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/graph_concepts.hpp"

// Forward-declare ImplicitBoundsGraph
namespace Scine {
namespace molassembler {
namespace distance_geometry {
class ImplicitBoundsGraph;
class ExplicitBoundsGraph;
} // namespace distance_geometry
} // namespace molassembler
} // namespace Scine

namespace boost {

template<
  typename VertexDescriptor,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap
>
void gor1_scan_helper(
  const VertexDescriptor& vertex,
  const VertexDescriptor& targetVertex,
  const double vertexDistance,
  const double edgeWeight,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map,
  std::stack<VertexDescriptor>& B
) {
  using ColorValue = typename property_traits<ColorMap>::value_type;
  using Color = color_traits<ColorValue>;

  if(vertexDistance + edgeWeight < get(distance_map, targetVertex)) {
    // Set distance for target vertex
    put(distance_map, targetVertex, vertexDistance + edgeWeight);

    // Mark target as vertex's descendant
    put(predecessor_map, targetVertex, vertex);

    // If not marked scanned, add it to B and mark it
    if(get(color_map, targetVertex) != Color::black()) {
      B.push(targetVertex);
      put(color_map, targetVertex, Color::gray());
    }
  }
}

template<
  typename VertexDescriptor,
  class IncidenceGraph,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap
>
std::enable_if_t<
  std::is_same<IncidenceGraph, Scine::molassembler::distance_geometry::ImplicitBoundsGraph>::value,
  void
> gor1_ig_scan(
  const VertexDescriptor& vertex,
  const IncidenceGraph& graph,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map,
  std::stack<VertexDescriptor>& B
) {
  auto vertexDistance = get(distance_map, vertex);

  /* ig edge_iterator needs to have local state of showing implicit edges or
   * not so that out_edges can be supplied with a boolean of show_implicit or
   * not
   *
   * if the (current!) shortest paths distance to the vertex is greater than
   * the largest implicit vdw bound, supply false
   *
   * the idea being that if we are messing around in the left graph and our
   * current shortest paths distance is higher than the largest implicit bound,
   * then we should have taken an implicit bound instead of an explicit in-group
   * arc instead of the last arc in the path
   *
   * Maybe just: If shortest paths distance > 0, disregard implicit edges
   *
   * Try optimized iterators for the cases
   * - left and shortest paths distance == 0
   *   -> all edges
   * - left and shortest paths distance > max implicit lower bound of vertex
   *   -> only left-left edges, no implicit or explicit left-right edges
   * - right: only right-right edges
   */


  if(IncidenceGraph::isLeft(vertex) && vertexDistance <= graph.maximalImplicitLowerBound(vertex)) {
    auto out_iter_pair = out_edges(vertex, graph);

    while(out_iter_pair.first != out_iter_pair.second) {
      VertexDescriptor targetVertex = out_iter_pair.first.target();
      double edgeWeight = out_iter_pair.first.weight();

      gor1_scan_helper(
        vertex,
        targetVertex,
        vertexDistance,
        edgeWeight,
        predecessor_map,
        color_map,
        distance_map,
        B
      );

      ++out_iter_pair.first;
    }
  } else {
    auto iter = graph.in_group_edges_begin(vertex);
    const auto end = graph.in_group_edges_end(vertex);

    while(iter != end) {
      VertexDescriptor targetVertex = iter.target();
      double edgeWeight = iter.weight();

      gor1_scan_helper(
        vertex,
        targetVertex,
        vertexDistance,
        edgeWeight,
        predecessor_map,
        color_map,
        distance_map,
        B
      );

      ++iter;
    }
  }
}


template<
  class IncidenceGraph,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap,
  typename VertexDescriptor
>
std::enable_if_t<
  std::is_same<IncidenceGraph, Scine::molassembler::distance_geometry::ImplicitBoundsGraph>::value,
  bool
> gor1_ig_shortest_paths(
  const IncidenceGraph& graph,
  const VertexDescriptor& root_vertex,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map
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

  while(!B.empty()) {
    // Compute A from B, emptying B in the process
    while(!B.empty()) {
      VertexDescriptor v = B.top();
      B.pop();

      ColorValue v_color = get(color_map, v);

      if(v_color == Color::black()) {
        A.push(v);
      } else if(v_color == Color::gray()) {
        /* In order for this vertex to be between any parents and predecessors
         * in B's and hence also A's stack order, must re-push B before adding
         * any of v's descendants to B. Mark it black though to avoid re-doing
         * this when the vertex is popped from the stack a second time
         */
        B.push(v);
        put(color_map, v, Color::black());

        gor1_ig_scan(
          v,
          graph,
          predecessor_map,
          color_map,
          distance_map,
          B
        );
      }
    }

    // Scan all elements in A, populating B with nodes added to the tree
    while(!A.empty()) {
      VertexDescriptor v = A.top();
      A.pop();

      // Scan
      gor1_ig_scan(
        v,
        graph,
        predecessor_map,
        color_map,
        distance_map,
        B
      );

      // Mark white
      put(color_map, v, Color::white());
    }
  }

  return true;
}

template<
  typename VertexDescriptor,
  class GraphClass,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap
>
void gor1_eg_scan(
  const VertexDescriptor& vertex,
  const GraphClass& graphWrapper,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map,
  std::stack<VertexDescriptor>& B
) {
  using ColorValue = typename property_traits<ColorMap>::value_type;
  using Color = color_traits<ColorValue>;

  const auto& graph = graphWrapper.graph();

  auto vertexDistance = get(distance_map, vertex);

  /* ig edge_iterator needs to have local state of showing implicit edges or
   * not so that out_edges can be supplied with a boolean of show_implicit or
   * not
   *
   * if the (current!) shortest paths distance to the vertex is greater than
   * the largest implicit vdw bound, supply false
   *
   * the idea being that if we are messing around in the left graph and our
   * current shortest paths distance is higher than the largest implicit bound,
   * then we should have taken an implicit bound instead of an explicit in-group
   * arc instead of the last arc in the path
   *
   * Maybe just: If shortest paths distance > 0, disregard implicit edges
   *
   * Try optimized iterators for the cases
   * - left and shortest paths distance == 0
   *   -> all edges
   * - left and shortest paths distance > max implicit lower bound of vertex
   *   -> only left-left edges, no implicit or explicit left-right edges
   * - right: only right-right edges
   */

  bool left = GraphClass::isLeft(vertex);

  if(
    left
    && vertexDistance <= graphWrapper.maximalImplicitLowerBound(vertex)
  ) {
    auto out_iter_pair = out_edges(vertex, graph);

    while(out_iter_pair.first != out_iter_pair.second) {
      auto edgeDescriptor = *out_iter_pair.first;
      VertexDescriptor targetVertex = target(edgeDescriptor, graph);
      double edgeWeight = get(edge_weight, graph, edgeDescriptor);

      if(vertexDistance + edgeWeight < get(distance_map, targetVertex)) {
        // Set distance for target vertex
        put(distance_map, targetVertex, vertexDistance + edgeWeight);

        // Mark target as vertex's descendant
        put(predecessor_map, targetVertex, vertex);

        // If not marked scanned, add it to B and mark it
        if(get(color_map, targetVertex) != Color::black()) {
          B.push(targetVertex);
          put(color_map, targetVertex, Color::gray());
        }
      }

      ++out_iter_pair.first;
    }
  } else {
    auto out_iter_pair = out_edges(vertex, graph);

    while(out_iter_pair.first != out_iter_pair.second) {
      auto edgeDescriptor = *out_iter_pair.first;
      VertexDescriptor targetVertex = target(edgeDescriptor, graph);

      if(left == GraphClass::isLeft(targetVertex)) {
        double edgeWeight = get(edge_weight, graph, edgeDescriptor);

        if(vertexDistance + edgeWeight < get(distance_map, targetVertex)) {
          // Set distance for target vertex
          put(distance_map, targetVertex, vertexDistance + edgeWeight);

          // Mark target as vertex's descendant
          put(predecessor_map, targetVertex, vertex);

          // If not marked scanned, add it to B and mark it
          if(get(color_map, targetVertex) != Color::black()) {
            B.push(targetVertex);
            put(color_map, targetVertex, Color::gray());
          }
        }
      }

      ++out_iter_pair.first;
    }
  }
}
template<
  class GraphClass,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap,
  typename VertexDescriptor
>
std::enable_if_t<
  std::is_same<GraphClass, Scine::molassembler::distance_geometry::ExplicitBoundsGraph>::value,
  bool
> gor1_eg_shortest_paths(
  const GraphClass& graphWrapper,
  const VertexDescriptor& root_vertex,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map
) {
  using IncidenceGraph = typename GraphClass::Graph;
  const auto& graph = graphWrapper.graph();

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

  while(!B.empty()) {
    // Compute A from B, emptying B in the process
    while(!B.empty()) {
      VertexDescriptor v = B.top();
      B.pop();

      ColorValue v_color = get(color_map, v);

      if(v_color == Color::black()) {
        A.push(v);
      } else if(v_color == Color::gray()) {
        /* In order for this vertex to be between any parents and predecessors
         * in B's and hence also A's stack order, must re-push B before adding
         * any of v's descendants to B. Mark it black though to avoid re-doing
         * this when the vertex is popped from the stack a second time
         */
        B.push(v);
        put(color_map, v, Color::black());

        gor1_eg_scan(
          v,
          graphWrapper,
          predecessor_map,
          color_map,
          distance_map,
          B
        );
      }
    }

    // Scan all elements in A, populating B with nodes added to the tree
    while(!A.empty()) {
      VertexDescriptor v = A.top();
      A.pop();

      // Scan
      gor1_eg_scan(
        v,
        graphWrapper,
        predecessor_map,
        color_map,
        distance_map,
        B
      );

      // Mark white
      put(color_map, v, Color::white());
    }
  }

  return true;
}

} // namespace boost

#endif
