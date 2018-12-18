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
  Scine::Utils::PositionCollection pos,
  const LengthUnit lengthUnit
) : positions(std::move(pos)) {
  if(lengthUnit == LengthUnit::Bohr) {
    positions *= Scine::Utils::Constants::angstrom_per_bohr;
  }
}

Scine::Utils::PositionCollection AngstromWrapper::getBohr() {
  if(_invalidated) {
    throw std::logic_error("AngstromWrapper was invalidated!");
  }

  _invalidated = true;
  positions *= Scine::Utils::Constants::bohr_per_angstrom;

  return std::move(positions);
}

} // namespace molassembler

} // namespace Scine
