/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Generate non-superposable atom-centered stereopermutations
 *
 * Main entry point into the library. From here, a set of rotationally unique
 * stereopermutations for a specific stereopermutation can be generated, with or without
 * counting the relative weights of unique stereopermutations.
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_GENERATE_UNIQUES_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_GENERATE_UNIQUES_H

#include "stereopermutation/Stereopermutation.h"

namespace Scine {

namespace stereopermutation {

/*! @brief Whether a stereopermutation has trans arranged linked substituents
 *
 * @complexity{@math{O(L)}}
 */
bool hasTransArrangedPairs(
  const Stereopermutation& stereopermutation,
  Symmetry::Shape shape
);

/*!
 * @brief Generate the set of rotationally unique stereopermutations for a given
 *   stereopermutation.
 *
 * By default does not remove trans-spanning groups (where a linked group's
 * directly bonded atoms span an angle of 180°).
 *
 * @note Gives NO guarantees as to satisfiability (if stereopermutation can be
 * fulfilled with real ligands)
 *
 * E.g. M (A-A)_3 generates a trans-trans-trans stereopermutation, which is extremely
 * hard to find actual ligands for that work.
 *
 * The satisfiability of stereopermutation must be checked before trying to embed
 * structures with completely nonsensical constraints. Perhaps restrict A-A
 * ligands with bridge length 4 (chelating atoms included), maybe even up to 6
 * to cis arrangements. Xantphos (with bridge length 7) is the smallest
 * trans-spanning ligand mentioned in Wikipedia.
 *
 * @complexity{@math{O(S!)} where @math{S} is the size of the involved symmetry}
 */
std::vector<Stereopermutation> uniques(
  const Stereopermutation& initial,
  Symmetry::Shape shape,
  bool removeTransSpanningGroups = false
);

//! Data class for uniqueStereopermutations including weights
struct StereopermutationsWithWeights {
  //! Stereopermutations found
  std::vector<Stereopermutation> stereopermutations;
  //! Relative statistical occurence of the stereopermutations
  std::vector<unsigned> weights;
};

/*!
 * Returns the set of rotationally unique stereopermutations including absolute
 * occurrence counts.
 *
 * See uniques()
 */
StereopermutationsWithWeights uniquesWithWeights(
  const Stereopermutation& initial,
  Symmetry::Shape shape,
  bool removeTransSpanningGroups = false
);

} // namespace stereopermutation

} // namespace Scine

#endif
