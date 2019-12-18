/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/AngstromWrapper.h"

#include "Utils/Constants.h"

namespace Scine {

namespace molassembler {

AngstromWrapper::AngstromWrapper(const unsigned N)
  : positions(Utils::PositionCollection::Zero(N, 3)) {}

AngstromWrapper::AngstromWrapper(
  const Utils::PositionCollection& pos,
  const LengthUnit lengthUnit
) {
  if(lengthUnit == LengthUnit::Bohr) {
    positions = Utils::Constants::angstrom_per_bohr * pos;
  } else {
    positions = pos;
  }
}

Utils::PositionCollection AngstromWrapper::getBohr() const {
  return positions * Utils::Constants::bohr_per_angstrom;
}

} // namespace molassembler

} // namespace Scine
