#ifndef INCLUDE_MOLASSEMBLER_DESCRIPTORS_H
#define INCLUDE_MOLASSEMBLER_DESCRIPTORS_H

/*!@file
 *
 * Contains various descriptor calculations for Molecule instances
 */

namespace molassembler {

// Forward-declare molcule
class Molecule;

/*! Calculates number of freely rotatable bonds in a molecule
 *
 * Criteria:
 * - Must be a single bond
 * - Neither atom connected by the bond may be terminal
 * - May not be member of a small cycle. The threshold of this criterion
 *   can be set by the user. An edge that is a member of a cycle smaller than
 *   the threshold is not considered rotatable. The default is 5. To deactivate
 *   this feature, pass std::numeric_limits<unsigned>::max().
 *
 */
unsigned numRotatableBonds(
  const Molecule& mol,
  unsigned cycleThreshold = 5
);

} // namespace molassembler

#endif
