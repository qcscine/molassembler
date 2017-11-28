#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_utility.hpp"
#include "boost/graph/isomorphism.hpp"

using GraphType = boost::adjacency_list<
  /* OutEdgeListS = Type of Container for edges of a vertex
   * Options: vector, list, slist, set, multiset, unordered_set
   * Choice: setS, enforces absence of parallel edges in graph
   */
  boost::setS,
  /* VertexListS = Type of Container for vertices
   * Options: vector, list, slist, set, multiset, unordered_set
   * Choice: vecS, removing vertices is rare, keep memory use limited
   * Consequence: operation remove_vertex() invalidates:
   *   - Vertex descriptors / iterators
   *   - Edge descriptors / iterators
   *   - Adjacency iterators
   *
   *   Upshot is that graph traversal is faster
   */
  boost::vecS,
  /* DirectedS = Is the graph directed or not?
   * Choice: Undirected
   */
  boost::undirectedS
>;
