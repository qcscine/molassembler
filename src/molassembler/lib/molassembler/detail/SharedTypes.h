#ifndef INCLUDE_MOLASSEMBLER_SHARED_TYPES_H
#define INCLUDE_MOLASSEMBLER_SHARED_TYPES_H

// External libraries
#include "boost/graph/adjacency_list.hpp"

// In-house libraries
#include "Delib/ElementTypes.h"

#include <limits>

/*! @file
 *
 * Central types required across the entire project are defined here.
 */

namespace molassembler {

/* Global typedefs */
/*!
 * Bond type enumeration. Besides the classic organic single, double and triple
 * bonds, bond orders up to sextuple are explicitly included.
 *
 * Although currently unused, Aromatic and Eta bonds are included in
 * anticipation of their necessity.
 */
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

enum class LengthUnit {
  Bohr,
  Angstrom
};


//! Boost graph vertex and edge property types
namespace GraphDetail {

struct VertexData {
  Delib::ElementType elementType;
};

struct EdgeData {
  BondType bondType;
};

} // namespace GraphDetail

/*!
 * The type of the molecular graph. An adjacency list is used due to sparsity
 * of the molecular graph. Further explanation of graph representation choices
 * are found in the code.
 */
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

//! Shorthand to boost graph vertex descriptor
using AtomIndexType = GraphType::vertex_descriptor;

//! Shorthand to boost graph edge descriptor
using EdgeIndexType = GraphType::edge_descriptor;

//! Descriptive name for dlib indices
using dlibIndexType = long;

/*! For bitmasks grouping components of immediate atom environments
 *
 * Differing strictnesses of comparisons may be desirable for various
 * purposes, hence a modular comparison function is provided.
 */
enum class AtomEnvironmentComponents : unsigned {
  ElementTypes,
  BondOrders,
  Symmetries,
  Stereopermutations // Symmetries must be set in conjunction with this
};

namespace Stereocenters {

constexpr AtomIndexType removalPlaceholder = std::numeric_limits<AtomIndexType>::max();

} // namespace Stereocenters

} // namespace molassembler

#endif
