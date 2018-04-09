#include <stack>
#include <limits>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/graph_concepts.hpp"

namespace boost {

namespace detail {

struct DummyGor1Visitor {
/* Stack operations */
  template<typename VertexDescriptor, typename IncidenceGraph>
  void a_push(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  template<typename VertexDescriptor, typename IncidenceGraph>
  void a_pop(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  template<typename VertexDescriptor, typename IncidenceGraph>
  void b_push(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  template<typename VertexDescriptor, typename IncidenceGraph>
  void b_pop(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

/* Edge operations */
  template<typename EdgeDescriptor, typename IncidenceGraph>
  void examine_edge(const EdgeDescriptor& /* e */, const IncidenceGraph& /* g */) {}

  template<typename EdgeDescriptor, typename IncidenceGraph>
  void relax_edge(const EdgeDescriptor& /* e */, const IncidenceGraph& /* g */) {}

/* Vertex colorings */
  template<typename VertexDescriptor, typename IncidenceGraph>
  void mark_white(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  template<typename VertexDescriptor, typename IncidenceGraph>
  void mark_gray(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}

  template<typename VertexDescriptor, typename IncidenceGraph>
  void mark_black(const VertexDescriptor& /* v */, const IncidenceGraph& /* g */) {}
};

} // namespace detail

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

        gor1_simplified_scan(
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
      gor1_simplified_scan(
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
