/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Subgraph matching algorithms for molecules and graphs
 */

#ifndef INCLUDE_MOLASSEMBLER_SUBGRAPHS_H
#define INCLUDE_MOLASSEMBLER_SUBGRAPHS_H

#include "boost/bimap.hpp"
#include "Molassembler/Types.h"

namespace Scine {
namespace Molassembler {

// Forward-declarations
class Molecule;
class Graph;

namespace Subgraphs {

/*!
 * @brief Type used to represent mappings.
 *
 * Left vertices are needle atom indices, right vertices are haystack atom
 * indices
 */
using IndexMap = boost::bimap<AtomIndex, AtomIndex>;

/*! @brief Matching strictness for vertices
 *
 * Sets the conditions for which vertices can be considered matches.
 * Strictness increases with the value of the underlying enum and each
 * successive enum encompasses all previous comparison components.
 */
enum class MASM_EXPORT VertexStrictness : unsigned {
  //! Element type must be the same
  ElementType,
  /*!
   * There must be a low-effort transition from the shape of the needle
   * vertex to the larger shape of the haystack vertex.
   *
   * E.g. Trigonal pyramidal is subsumed in tetrahedral, seesaw is subsumed in
   * square pyramidal, square pyramidal is subsumed in octahedral, etc.
   *
   * @warning This strictness is not implemented
   */
  SubsumeShape,
  /*!
   * If a needle vertex carries an assigned stereopermutator, a haystack vertex
   * matches only if its chiral state encompasses
   *
   * @warning This strictness is not implemented
   */
  SubsumeStereopermutation
};

/**
 * @brief Matching strictness for edges
 *
 * Sets the conditions for which edges can be considered matches.
 * Strictness increases with the value of the underlying enum and each
 * successive enum encompasses all previous comparison components.
 */
enum class MASM_EXPORT EdgeStrictness : unsigned {
  //! No constraints are set upon edge matching besides graph topography
  Topographic,
  //! Bond types must match exactly
  BondType,
  /*!
   * If a needle bond carries an assigned stereopermutator, a haystack bond
   * must subsume its chiral state.
   *
   * E.g. Z-diazene (one assignment of two possible) in pyridazene (one
   * assignment of one possible), but not E-diazene in pyridazene (chiral state
   * is mismatched).
   *
   * @warning This strictness is not implemented
   */
  SubsumeStereopermutation
};

/**
 * @brief Searches for subgraphs of needle in haystack
 *
 * @param needle The smaller graph to search for
 * @param haystack The larger graph to search in
 * @params vertexStrictness Strictness with which to allow vertex matching.
 *   Maximum strictness for subgraph isomorphism is
 *   VertexStrictness::ElementType.
 * @params edgeStrictness Strictness with which to allow edge matching. Maximum
 *   strictness for subgraph isomorphism is EdgeStrictness::BondType.
 *
 * @return List of index mappings of vertices the needle to vertices of the
 *   haystack
 */
MASM_EXPORT std::vector<IndexMap> complete(
  const Graph& needle,
  const Graph& haystack,
  VertexStrictness vertexStrictness = VertexStrictness::ElementType,
  EdgeStrictness edgeStrictness = EdgeStrictness::Topographic
);

/**
 * @brief Searches for subgraphs of needle in haystack
 *
 * @param needle The smaller molecule to search for
 * @param haystack The larger molecule to search in
 * @params vertexStrictness Strictness with which to allow vertex matching.
 *   Maximum implemented strictness for molecule subgraph isomorphism is
 *   VertexStrictness::ElementType.
 * @params edgeStrictness Strictness with which to allow edge matching. Maximum
 *   implemented strictness for molecule subgraph isomorphism is
 *   EdgeStrictness::BondType.
 *
 * @return List of index mappings of vertices the needle to vertices of the
 *   haystack
 */
MASM_EXPORT std::vector<IndexMap> complete(
  const Molecule& needle,
  const Molecule& haystack,
  VertexStrictness vertexStrictness = VertexStrictness::ElementType,
  EdgeStrictness edgeStrictness = EdgeStrictness::Topographic
);

/*!
 * @brief Find mappings for the maximum common subgraph between two graphs
 *
 * Finds an index mapping from a to b representing all found maximum common
 * subgraphs (if present).
 *
 * @params a The first graph
 * @params b The second graph
 * @params vertexStrictness Strictness with which to allow vertex matching.
 *   Maximum strictness for graph MCS is VertexStrictness::ElementType.
 * @params edgeStrictness Strictness with which to allow edge matching. Maximum
 *   strictness for graph MCS is EdgeStrictness::BondType.
 *
 * @complexity{@math{O(N_1 \cdot N_2)} where @math{N_i} is the number of
 * vertices in graph @math{i}}
 *
 * @warning For subgraph comparison, only element and bond types are considered.
 * Stereocenters and Stereopermutations are not graph-local properties suitable
 * to a common substructure matching.
 */
MASM_EXPORT std::vector<IndexMap> maximum(
  const Graph& a,
  const Graph& b,
  VertexStrictness vertexStrictness = VertexStrictness::ElementType,
  EdgeStrictness edgeStrictness = EdgeStrictness::Topographic
);

/*!
 * @brief Find mappings for the maximum common subgraph between two molecules
 *
 * Finds an index mapping from a to b representing all found maximum common
 * subgraphs (if present).
 *
 * @params a The first molecule
 * @params b The second molecule
 * @params vertexStrictness Strictness with which to allow vertex matching.
 *   Maximum implemented strictness for molecule MCS is
 *   VertexStrictness::ElementType.
 * @params edgeStrictness Strictness with which to allow edge matching. Maximum
 *   implemented strictness for molecule MCS is EdgeStrictness::BondType.
 *
 * @complexity{@math{O(N_1 \cdot N_2)} where @math{N_i} is the number of
 * vertices in graph @math{i}}
 *
 * @warning For subgraph comparison, only element and bond types are considered.
 * Stereocenters and Stereopermutations are not graph-local properties suitable
 * to a common substructure matching.
 */
MASM_EXPORT std::vector<IndexMap> maximum(
  const Molecule& a,
  const Molecule& b,
  VertexStrictness vertexStrictness = VertexStrictness::ElementType,
  EdgeStrictness edgeStrictness = EdgeStrictness::Topographic
);

} // namespace Subgraphs
} // namespace Molassembler
} // namespace Scine

#endif
