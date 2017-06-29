#ifndef INCLUDE_COMMON_TYPEDEFS_H
#define INCLUDE_COMMON_TYPEDEFS_H

// External libraries
#include "boost/graph/adjacency_list.hpp"

// In-house libraries
#include "Delib/ElementTypes.h"

namespace MoleculeManip {

/* Global typedefs */
enum class BondType : unsigned {
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Sextuple,
    Aromatic,
    Eta
};

namespace GraphDetail {
  struct VertexData {
    Delib::ElementType elementType;
  };

  struct EdgeData {
    BondType bondType;
  };
} // namespace GraphDetail

//! The type of the molecular graph
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
  boost::undirectedS,
  /* VertexProperty = What information is stored about vertices?
   * Choice: Atom, containing an index and an element type
   */
  GraphDetail::VertexData,
  /* EdgeProperty = What information is stored about edges?
   * Choice: BondType, a custom enum class 
   */
  GraphDetail::EdgeData
  /* GraphProperty
   * Omitted, defaults
   */
  /* EdgeListS
   * Ommitted, defaults
   */
>;

using AtomIndexType = GraphType::vertex_descriptor;
using EdgeIndexType = GraphType::edge_descriptor;

using dlibIndexType = long;

} // namespace MoleculeManip

#endif
