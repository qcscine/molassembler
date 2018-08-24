#ifndef INCLUDE_MOLASSEMBLER_INNER_GRAPH_H
#define INCLUDE_MOLASSEMBLER_INNER_GRAPH_H

#include "boost/graph/adjacency_list.hpp"
#include "Delib/ElementTypes.h"

#include "molassembler/Layers/Types.h"

#include <limits>

namespace molassembler {

class InnerGraph {
public:
//!@name Member types
//!@{
  //! Information stored at each graph vertex
  struct VertexData {
    Delib::ElementType elementType;
  };

  //! Information stored at each graph edge
  struct EdgeData {
    BondType bondType;
  };

  /*! The type of the BGL graph.
   *
   * An adjacency list is used due to sparsity of the molecular graph.
   */
  using BGLType = boost::adjacency_list<
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
     */
    boost::vecS,
    // Molecular graphs are undirected
    boost::undirectedS,
    // VertexProperty: What information is stored about vertices?
    VertexData,
    // EdgeProperty: What information is stored about edges?
    EdgeData
    // GraphProperty: Omitted, defaults
    // EdgeListS: Omitted, defaults
  >;

  using Vertex = BGLType::vertex_descriptor;
  using Edge = BGLType::edge_descriptor;

  template<typename Iter>
  using Range = std::pair<Iter, Iter>;

  using VertexRange = Range<BGLType::vertex_iterator>;
  using EdgeRange = Range<BGLType::edge_iterator>;
  using AdjacentVertexRange = Range<BGLType::adjacency_iterator>;
  using IncidentEdgeRange = Range<BGLType::out_edge_iterator>;
//!@}

//!@name Static members
//!@{
  static constexpr Vertex removalPlaceholder = std::numeric_limits<Vertex>::max();
//!@}

//!@name Modification
//!@{
  BGLType& bgl();
//!@}

//!@name Information
//!@{
  Delib::ElementType elementType(Vertex a) const;
  BondType bondType(const Edge& edge) const;
  const BGLType& bgl() const;


  //! Precondition: Edge exists!
  Edge edge(Vertex a, Vertex b) const;
  Vertex source(const Edge& edge) const;
  Vertex target(const Edge& edge) const;
//!@}

//!@name Ranges
//!@{
  VertexRange vertices() const;
  EdgeRange edges() const;
  AdjacentVertexRange adjacents(Vertex a) const;
  IncidentEdgeRange edges(Vertex a) const;
//!@}

private:
  BGLType _graph;
};

} // namespace molassembler

#endif
