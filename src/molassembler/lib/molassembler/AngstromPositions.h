/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Wrapper class to typify Angstrom scale positional information
 *
 * Provides a class that exists purely to strongly separate position collections
 * in bohr units and angstrom units in library interfaces
 */

#ifndef INCLUDE_MOLASSEMBLER_ANGSTROM_POSITIONS_H
#define INCLUDE_MOLASSEMBLER_ANGSTROM_POSITIONS_H

#include "Utils/Typenames.h"

#include "molassembler/Types.h"

namespace Scine {
namespace molassembler {

/**
 * @brief A wrapper class around Utils' PositionCollection to emphasize
 *   that the positions stored therein are in Angstrom
 */
class MASM_EXPORT AngstromPositions {
public:
  //! Positions in angstrom units
  Utils::PositionCollection positions;

  //! Default ctor
  AngstromPositions() = default;
  //! Preallocate space for N positions
  explicit AngstromPositions(unsigned N);
  //! Convert from a Utils::PositionCollection
  explicit AngstromPositions(
    const Utils::PositionCollection& pos,
    LengthUnit lengthUnit = LengthUnit::Bohr
  );

  //! Fetch a bohr representation of the wrapped positions
  Utils::PositionCollection getBohr() const;
};

} // namespace molassmbler
} // namespace Scine

#endif
