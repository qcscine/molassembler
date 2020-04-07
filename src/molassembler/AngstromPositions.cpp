/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/AngstromPositions.h"

#include "Utils/Constants.h"

namespace Scine {
namespace molassembler {

AngstromPositions::AngstromPositions(const unsigned N)
  : positions(Utils::PositionCollection::Zero(N, 3)) {}

AngstromPositions::AngstromPositions(
  const Utils::PositionCollection& pos,
  const LengthUnit lengthUnit
) {
  if(lengthUnit == LengthUnit::Bohr) {
    positions = Utils::Constants::angstrom_per_bohr * pos;
  } else {
    positions = pos;
  }
}

Utils::PositionCollection AngstromPositions::getBohr() const {
  return positions * Utils::Constants::bohr_per_angstrom;
}

} // namespace molassembler
} // namespace Scine
