/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/AngstromWrapper.h"

#include "Utils/Constants.h"

namespace Scine {

namespace molassembler {

AngstromWrapper::AngstromWrapper(const unsigned N)
  : positions(Scine::Utils::PositionCollection::Zero(N, 3)) {}

AngstromWrapper::AngstromWrapper(
  const Scine::Utils::PositionCollection& pos,
  const LengthUnit lengthUnit
) {
  if(lengthUnit == LengthUnit::Bohr) {
    positions = Scine::Utils::Constants::angstrom_per_bohr * pos;
  } else {
    positions = pos;
  }
}

Scine::Utils::PositionCollection AngstromWrapper::getBohr() const {
  return positions * Scine::Utils::Constants::bohr_per_angstrom;
}

} // namespace molassembler

} // namespace Scine
