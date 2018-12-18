/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/AngstromWrapper.h"

#include "Delib/Constants.h"

namespace Scine {

namespace molassembler {

AngstromWrapper::AngstromWrapper(const unsigned N) : positions(N) {}

AngstromWrapper::AngstromWrapper(
  Delib::PositionCollection pos,
  const LengthUnit lengthUnit
) : positions {std::move(pos)} {
  if(lengthUnit == LengthUnit::Bohr) {
    for(auto& position : positions) {
      position *= Delib::angstrom_per_bohr;
    }
  }
}

Delib::PositionCollection AngstromWrapper::getBohr() {
  if(_invalidated) {
    throw std::logic_error("AngstromWrapper was invalidated!");
  }

  _invalidated = true;
  for(auto& position : positions) {
    position *= Delib::bohr_per_angstrom;
  }

  return positions;
}

} // namespace molassembler

} // namespace Scine
