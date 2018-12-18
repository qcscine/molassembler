/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Predicates to determine what kinds of stereoisomers molecule pairs are
 */

#ifndef INCLUDE_MOLASSEMBLER_ISOMERS_H
#define INCLUDE_MOLASSEMBLER_ISOMERS_H

namespace Scine {

namespace molassembler {

// Forward-declarations
class Molecule;

/**
 * @brief Tag for function dispatch indicating that two molecules are enumerated
 *   the same
 */
struct SameIndexingTag {};

/*!
 * @brief Returns whether three dimensional representations of two molecules
 *   are mirror images of one another
 *
 * This algorithm performs a graph isomorphism and tries to match atom
 * stereopermutators with its mirror assignments.
 *
 * @throws std::runtime_error Since it's unimplemented so far.
 * @todo implement
 */
bool enantiomeric(const Molecule& /* a */, const Molecule& /* b */);

/*!
 * @brief Returns whether three dimensional representations of two molecules
 *   are mirror images of one another if both molecules are indexed identically
 *
 * This algorithm looks through both Molecules' StereopermutatorLists and
 * checks whether each pair of assigned AtomStereopermutators are mirror images
 * of one another.
 *
 * @pre The two molecules are indexed the same, i.e. a graph isomorphism
 *   based on element types, bond types and symmetries would yield an identity
 *   mapping. This is always the case if a molecule is created by copying and
 *   re-assigning stereopermutators.
 *
 * @note When determining enantiomerism of a pair of AtomStereopermutators, if
 *   either permutator is unassigned, then three-dimensional representations @b
 *   may be enantiomeric or not (determined by random choice of assignments at
 *   conformer generation time). In these cases, this function returns false.
 *
 * @throws std::logic_error If the molecules are not of the same size or the
 *   molecules are not indexed identically.
 *
 * @param a The first molecule to compare
 * @param b The second molecule to compare
 * @param sameIndexingTag Tag to indicate that the two Molecules are indexed
 *   identically
 *
 * @returns Whether two Molecules' conformers are enantiomeric
 */
bool enantiomeric(
  const Molecule& a,
  const Molecule& b,
  SameIndexingTag /* sameIndexingTag */
);

/**
 * @brief Generates a molecule's enantiomer
 *
 * @param a The molecule whose enantiomer to generate
 *
 * @note If there are no atom stereopermutators in this molecule with more than
 *   one stereopermutations, yields a Molecule identical to @p a
 *
 * @todo Untested
 *
 * @return The enantiomer to a molecule
 */
Molecule enantiomer(const Molecule& a);

} // namespace molassembler

} // namespace Scine

#endif
