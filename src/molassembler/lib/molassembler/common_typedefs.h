#ifndef INCLUDE_COMMON_TYPEDEFS_H
#define INCLUDE_COMMON_TYPEDEFS_H

// External libraries
#include "boost/graph/adjacency_list.hpp"

// In-house libraries
#include "Delib/ElementTypes.h"

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

/*!
 * Specifies for which temperature regime the Molecule is being modeled.
 * This currently affects only whether nitrogen atoms with a trigonal
 * pyramidal geometry are considered stereocenters. At low temperatures,
 * there is no pyramidal inversion, and these nitrogens are considered
 * valid stereocenters. At high temperatures, only nitrogen geometries
 * that are in particularly strained cycles (3, 4) are considered
 * stereocenters.
 */
enum class TemperatureRegime {Low, High};

/*!
 * Specifies the effects of graph modifications. In case a substituent is
 * added or removed at a stereocenter, an attempt is made to transfer chiral
 * information into the new geometry. How this attempt is made can be altered:
 * - None: Don't try at all (complete loss of information)
 * - Effortless and unique: Use only completely unambiguous zero-effort
 *   mappings. Those are the green edges in the graphs below. Note that the
 *   ligand gain situation from square planar to square pyramidal is not
 *   unique, and therefore not shown as green.
 * - Unique: Propagates if the best symmetry mapping is unique, i.e. there
 *   are no other mappings with the same quality measures. This enables all
 *   green and black edges.
 * - RandomFromMultipleBest: Chooses randomly from the set of best mappings,
 *   permitting chiral state propagation in all cases. So propagating chiral
 *   state from square planar to square pyramidal is now possible -- there
 *   are two ways of placing the new apical ligand -- but you only get one of
 *   them.
 *
 * The following graphs are to illustrate which chiral information transfers
 * are possible under which setting.
 *
 * The first is for the situation of ligand gain. There are very few edges and
 * symmetries in this graph because in nearly all cases, there are multiple
 * best mappings, and propagation in these cases are enabled only in
 * RandomFromMultipleBest, in which case the graph gets quite dense, so it's
 * not shown.
 *
 * \dot
 *  digraph g {
 *    graph [fontname = "Arial", nodesep="0.2", ranksep="0.8"];
 *    node [fontname = "Arial", style = "filled", fillcolor="white"];
 *    edge [fontname = "Arial", penwidth=2, labelfontsize="10"];
 *    rankdir="LR";
 *
 *    trigonalpyramidal [label="trigonal pyramidal"];
 *    Tshaped [label="T-shaped"];
 *    tetrahedral [label="tetrahedral"];
 *    squareplanar [label="square planar"];
 *    seesaw [label="seesaw"];
 *    squarepyramidal [label="square pyramidal"];
 *    trigonalbipyramidal [label="trigonal bipyramidal"];
 *    octahedral [label="octahedral"];
 *    pentagonalpyramidal [label="pentagonal pyramidal"];
 *    pentagonalbipyramidal [label="pentagonal bipyramidal"];
 *    trigonalpyramidal -> tetrahedral [color="forestgreen"];
 *    Tshaped -> squareplanar [color="forestgreen"];
 *    seesaw -> trigonalbipyramidal [color="forestgreen"];
 *    squarepyramidal -> octahedral [color="forestgreen"];
 *    pentagonalpyramidal -> pentagonalbipyramidal [color="forestgreen"];
 *  }
 * \enddot
 *
 * The second is for the situation of ligand loss. The edge label indicates
 * the group of symmetry ligands from which a ligand is being removed.
 *
 * \dot
 *  digraph g {
 *    graph [fontname = "Arial", nodesep="0.8", ranksep="0.8"];
 *    node [fontname = "Arial", style = "filled", fillcolor="white"];
 *    edge [fontname = "Arial", penwidth=2, labelfontsize="10"];
 *    rankdir="LR";
 *
 *    trigonalpyramidal [label="trigonal pyramidal"];
 *    Tshaped [label="T-shaped"];
 *    tetrahedral [label="tetrahedral"];
 *    squareplanar [label="square planar"];
 *    seesaw [label="seesaw"];
 *    squarepyramidal [label="square pyramidal"];
 *    trigonalbipyramidal [label="trigonal bipyramidal"];
 *    pentagonalplanar [label="pentagonal planar"];
 *    octahedral [label="octahedral"];
 *    trigonalprismatic [label="trigonal prismatic"];
 *    pentagonalpyramidal [label="pentagonal pyramidal"];
 *    pentagonalbipyramidal [label="pentagonal bipyramidal"];
 *
 *    tetrahedral -> trigonalpyramidal [color="forestgreen", label="any"];
 *    squareplanar -> Tshaped [color="forestgreen", label="any"];
 *    seesaw -> trigonalpyramidal [color="black", label="axial"];
 *    seesaw -> Tshaped [color="black", label="1 axial"];
 *    seesaw -> Tshaped [color="forestgreen", label="equat."];
 *    squarepyramidal -> tetrahedral [color="black", label="equat."];
 *    squarepyramidal -> squareplanar [color="black", label="equat."];
 *    squarepyramidal -> squareplanar [color="forestgreen", label="axial"];
 *    squarepyramidal -> seesaw [color="black", label="equat."];
 *    trigonalbipyramidal -> tetrahedral [color="black", label="equat."];
 *    trigonalbipyramidal -> tetrahedral [color="black", label="axial"];
 *    trigonalbipyramidal -> squareplanar [color="black", label="equat."];
 *    trigonalbipyramidal -> seesaw [color="forestgreen", label="equat."];
 *    pentagonalplanar -> squareplanar [color="black", label="any"];
 *    octahedral -> squarepyramidal [color="forestgreen", label="any"];
 *    trigonalprismatic -> pentagonalplanar [color="black", label="any"];
 *    pentagonalpyramidal -> squarepyramidal [color="black", label="equat."];
 *    pentagonalpyramidal -> pentagonalplanar [color="black", label="equat."];
 *    pentagonalpyramidal -> pentagonalplanar [color="forestgreen", label="apical"];
 *    pentagonalbipyramidal -> octahedral [color="black", label="equat."];
 *    pentagonalbipyramidal -> pentagonalpyramidal [color="forestgreen", label="axial"];
 *  }
 *   \enddot
 */
enum class ChiralStatePreservation {
  None,
  EffortlessAndUnique,
  Unique,
  RandomFromMultipleBest
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

} // namespace molassembler

#endif
