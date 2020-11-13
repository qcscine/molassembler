/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Calculates properties of molecules
 */

#ifndef INCLUDE_MOLASSEMBLER_DESCRIPTORS_H
#define INCLUDE_MOLASSEMBLER_DESCRIPTORS_H

#include "Molassembler/Types.h"
#include <vector>

namespace Scine {
namespace Molassembler {

// Forward-declare molcule
class Molecule;

/*! @brief Calculates a number of freely rotatable bonds in a molecule
 *
 * The number of rotatable bonds is calculated as a bond-wise sum of
 * contributions that is rounded at the end:
 *
 * A bond can contribute to the number of rotatable bonds if
 * - It is of bond type Single
 * - Neither atom connected by the bond is terminal
 * - There is no assigned BondStereopermutator on the bond
 *
 * A bond meeting the prior criteria:
 * - If not part of a cycle, contributes a full rotatable bond to the sum
 * - If part of a cycle, contributes (S - 3) / S to the sum (where S is
 *   the cycle size).
 *
 * @complexity{@math{\Theta(B)} where @math{B} is the number of bonds}
 *
 * @warning The number of rotatable bonds is an unphysical descriptor and
 *   definitions differ across libraries. Take the time to read the algorithm
 *   description implemented here and do some testing. If need be, all
 *   information used by this algorithm is accessible from the Molecule
 *   interface, and a custom algorithm can be implemented.
 */
MASM_EXPORT unsigned numRotatableBonds(const Molecule& mol);

/**
 * @brief Yields group memberships for ranking equivalent atoms
 *
 * Uses atom-stereopermutator ranking information to determine which parts of
 * molecules are completely ranking-equivalent. Note that ranking symmetry at
 * bonds is not found by this method, yielding more distinct atoms than there
 * truly are in many cases, such as in ethane, cyclobutane or benzene.
 *
 * @return A list of length equal to the number of atoms in @p mol of group
 * indices. Atoms with the same group index are ranking equivalent.
 */
MASM_EXPORT std::vector<unsigned> rankingEquivalentGroups(const Molecule& mol);

/**
 * @brief Determines non-ranking equivalent atoms in a molecule
 *
 * Uses atom-stereopermutator ranking information to determine which parts of
 * molecules are completely ranking-equivalent. Note that ranking symmetry at
 * bonds is not found by this method, yielding more distinct atoms than there
 * truly are in many cases, such as in ethane, cyclobutane or benzene.
 *
 * @return Unordered list of non-ranking equivalent atoms
 */
MASM_EXPORT std::vector<AtomIndex> rankingDistinctAtoms(const Molecule& mol);

} // namespace Molassembler
} // namespace Scine

#endif
