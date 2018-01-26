#include <stack>
#include <limits>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/graph_concepts.hpp"

// Forward-declare ImplicitGraph
namespace MoleculeManip {
namespace DistanceGeometry {
class ImplicitGraph;
}
}

namespace boost {

/* SPG edge_iterator needs to have local state of showing implicit edges or
 * not so that out_edges can be supplied with a boolean of show_implicit or
 * not
 *
 * if the (current!) shortest paths distance to the vertex is greater than
 * the largest implicit vdw bound, supply false
 *
 * the idea being that if we are fucking around in the left graph and our
 * current shortest paths distance is higher than the largest implicit bound,
 * then we should have taken an implicit bound instead of an explicit in-group
 * arc instead of the last arc in the path
 *
 * Maybe just: If shortest paths distance > 0, disregard implicit edges
 *
 * TODO Try optimized iterators for the cases
 * - left and shortest paths distance == 0
 *   -> all edges
 * - left and shortest paths distance > 0
 *   -> only left-left edges, no implicit or explicit left-right edges
 * - right: only right-right edges
 */

template<
  typename VertexDescriptor,
  class IncidenceGraph,
  class DistanceMap,
  class PredecessorMap,
  class ColorMap
> 
void gor1_simplified_scan(
  const VertexDescriptor& vertex,
  const IncidenceGraph& graph,
  PredecessorMap& predecessor_map,
  ColorMap& color_map,
  DistanceMap& distance_map,
  std::stack<VertexDescriptor>& B
) {
  using ColorValue = typename property_traits<ColorMap>::value_type;
  using Color = color_traits<ColorValue>;

  auto vertexDistance = get(distance_map, vertex);
  auto out_iter_pair = out_edges(vertex, graph);

  while(out_iter_pair.first != out_iter_pair.second) {
    VertexDescriptor targetVertex = target(*out_iter_pair.first, graph);

    double edgeWeight = get(boost::edge_weight, graph, *out_iter_pair.first);

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

        gor1_simplified_scan(
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
      gor1_simplified_scan(
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

} // namespace boost
