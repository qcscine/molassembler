#ifndef INCLUDE_MOLASSEMBLER_DESCRIPTORS_H
#define INCLUDE_MOLASSEMBLER_DESCRIPTORS_H

/*!@file
 *
 * Contains various descriptor calculations for Molecule instances
 */

namespace molassembler {

// Forward-declare molcule
class Molecule;

/*! Calculates a number of freely rotatable bonds in a molecule
 *
 * A freely rotatable bond:
 * - Must be a single bond
 * - Neither atom connected by the bond may be terminal
 * - There may not be an assigned stereogenic BondStereocenter on the bond
 *
 * A rotatable bond in a cycle:
 * - If the edge is part of a cycle, it contributes at most (S - 3) / 3 to the
 *   total number of rotatable bonds. The final (fractional) count is rounded.
 */
unsigned numRotatableBonds(const Molecule& mol);

} // namespace molassembler

#endif
