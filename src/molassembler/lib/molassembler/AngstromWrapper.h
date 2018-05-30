#ifndef INCLUDE_MOLASSEMBLER_ANGSTROM_POSITIONS_H
#define INCLUDE_MOLASSEMBLER_ANGSTROM_POSITIONS_H

#include "Delib/PositionCollection.h"
#include "Delib/Constants.h"
#include "common_typedefs.h"

namespace molassembler {

class AngstromWrapper {
private:
  bool invalidated = false;

public:
  Delib::PositionCollection positions;

  AngstromWrapper() = default;
  inline explicit AngstromWrapper(const unsigned N) : positions(N) {}
  inline explicit AngstromWrapper(
    Delib::PositionCollection pos,
    const LengthUnit lengthUnit = LengthUnit::Bohr
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
  inline Delib::PositionCollection getBohr() {
    if(invalidated) {
      throw std::logic_error("AngstromWrapper was invalidated!");
    }

    invalidated = true;
    for(auto& position : positions) {
      position *= Delib::bohr_per_angstrom;
    }

    return positions;
  }
};

} // namespace molassmbler

#endif
