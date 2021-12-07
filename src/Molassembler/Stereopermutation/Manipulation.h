/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Manipulations of stereopermutations with shapes
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_MANIPULATION_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_MANIPULATION_H

#include "Molassembler/Stereopermutation/Stereopermutation.h"
#include <unordered_set>

namespace Scine {
namespace Molassembler {
namespace Stereopermutations {

using UnorderedStereopermutations = std::unordered_set<Stereopermutation, boost::hash<Stereopermutation>>;

/*! @brief Generate all superimposable rotations of a stereopermutation
 *
 * Generates a set of all rotational equivalents of this Stereopermutation as
 * defined by its shape template parameter.
 *
 * @complexity{@math{O(\prod_i^Rm_i)} where @math{R} is the set of rotations and
 * @math{m_i} is the multiplicity of rotation @math{i}}
 */
MASM_EXPORT std::vector<Stereopermutation> generateAllRotations(
  Stereopermutation s,
  Shapes::Shape shape
);

/*! @brief whether this Stereopermutation is rotationally superimposable with another.
 *
 * @complexity{As generateAllRotations}
 */
MASM_EXPORT bool rotationallySuperimposable(
  Stereopermutation a,
  const Stereopermutation& b,
  Shapes::Shape shape
);

/*!
 * @brief Checks whether a stereopermutation is a mirror image of another
 *   within a particular shape
 *
 * @complexity{As generateAllRotations}
 *
 * @returns boost::none If the shape does not generate enantiomers
 * @returns true If the shape has enantiomers, and this is the enantiomeric
 *   stereopermutation to @p other
 * @return false If the shape has enantiomers, and this is not the
 *   enantiomeric to @p other
 */
MASM_EXPORT boost::optional<bool> enantiomer(
  const Stereopermutation& a,
  const Stereopermutation& b,
  Shapes::Shape shape
);

/*! @brief Whether a stereopermutation has trans arranged linked substituents
 *
 * @complexity{@math{O(L)}}
 */
MASM_EXPORT bool hasTransArrangedLinks(
  const Stereopermutation& s,
  Shapes::Shape shape
);

struct MASM_EXPORT Uniques {
  std::vector<Stereopermutation> list;
  std::vector<unsigned> weights;
};

/*!
 * @brief Generate the set of rotationally unique stereopermutations for a
 *   given stereopermutation.
 *
 * By default does not remove trans-spanning groups (where a linked group's
 * directly bonded atoms span an angle of 180°).
 *
 * @note Gives NO guarantees as to satisfiability (if stereopermutation can be
 * fulfilled with real ligands)
 *
 * E.g. M (A-A)_3 generates a trans-trans-trans stereopermutation, which is
 * extremely hard to find actual ligands for that work.
 *
 * The satisfiability of a stereopermutation must be checked before trying to
 * embed structures with completely nonsensical constraints.
 *
 * @complexity{@math{O(S!)} where @math{S} is the size of the involved symmetry}
 */
MASM_EXPORT Uniques uniques(
  const Stereopermutation& base,
  Shapes::Shape shape,
  bool removeTransSpanningGroups = false
);

} // namespace Stereopermutations
} // namespace Molassembler
} // namespace Scine

#endif
