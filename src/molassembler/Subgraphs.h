/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Maximum common subgraph matching algorithms for molecules
 */

#ifndef INCLUDE_MOLASSEMBLER_SUBGRAPHS_H
#define INCLUDE_MOLASSEMBLER_SUBGRAPHS_H

#include "boost/bimap.hpp"
#include "molassembler/Types.h"

/*!@file
 *
 * Enables common subgraph computations
 */

namespace Scine {
namespace molassembler {

// Forward-declare Molecule
class Molecule;

namespace subgraphs {

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
   */
  SubsumeShape,
  /*!
   * If a needle vertex carries an assigned stereopermutator, a haystack vertex
   * matches only if its chiral state encompasses
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
  //! No constraints are set upon vertex matching besides graph topography
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
   */
  SubsumeStereopermutation
};

/*!
 * @brief Find mappings for the maximum common subgraph between two molecules
 *
 * Finds an index mapping from a to b representing all found maximum common
 * subgraphs (if present).
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
  EdgeStrictness edgeStrictness = EdgeStrictness::Topographic,
  bool removeHydrogenPermutations = true
);

} // namespace subgraphs
} // namespace molassembler
} // namespace Scine

#endif
