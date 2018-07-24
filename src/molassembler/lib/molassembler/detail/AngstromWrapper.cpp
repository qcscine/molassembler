#include "detail/AngstromWrapper.h"

#include "Delib/Constants.h"

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

/* TODO Unsure if I should explicitly move from this, or whether that would
 * be a pessimization in all cases, but it most clearly states intent - after
 * calling this function, AngstromWrapper is not supposed to be reused.
 */
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
