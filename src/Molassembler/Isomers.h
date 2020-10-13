/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Predicates to determine what kinds of stereoisomers molecule pairs are
 */

#ifndef INCLUDE_MOLASSEMBLER_ISOMERS_H
#define INCLUDE_MOLASSEMBLER_ISOMERS_H

#include "Molassembler/Export.h"

namespace Scine {
namespace Molassembler {

// Forward-declarations
class Molecule;

/*!
 * @brief Returns whether three dimensional representations of two molecules
 *   are mirror images of one another
 *
 * This algorithm looks through both Molecules' StereopermutatorLists and
 * checks whether each pair of assigned AtomStereopermutators are mirror images
 * of one another and there is at least one enantiomeric pair.
 *
 * @note When determining enantiomerism of a pair of AtomStereopermutators, if
 *   either permutator is unassigned, then three-dimensional representations @b
 *   may be enantiomeric or not (determined by random choice of assignments at
 *   conformer generation time). In these cases, this function returns false.
 */
MASM_EXPORT bool enantiomeric(
  const Molecule& a,
  const Molecule& b
);

/**
 * @brief Determine whether molecules are diastereomers of one another
 *
 * Two molecules are diastereomers if they are not mirror images of one another
 * and have different configurations at one or more stereocenters.
 */
MASM_EXPORT bool diastereomeric(
  const Molecule& a,
  const Molecule& b
);

/**
 * @brief Determine whether molecules are epimers of one another
 *
 * Two molecules are epimers if they differ in exactly one stereocenter.
 */
MASM_EXPORT bool epimeric(
  const Molecule& a,
  const Molecule& b
);

/** @brief Generates a molecule's enantiomer
 *
 * @complexity{@math{O(A)} where @math{A} is the number of chiral atom
 * stereopermutators of one molecule}
 *
 * @param a The molecule whose enantiomer to generate
 *
 * @note If there are no atom stereopermutators in this molecule with more than
 *   one stereopermutations, yields a Molecule identical to @p a
 *
 * @warning This has been tested very little
 *
 * @return The enantiomer to a molecule
 */
MASM_EXPORT Molecule enantiomer(const Molecule& a);

} // namespace Molassembler
} // namespace Scine

#endif
