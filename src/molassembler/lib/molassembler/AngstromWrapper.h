/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
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
 * @brief A wrapper class around Scine::Utils' PositionCollection to emphasize
 *   that the positions stored therein are in Angstrom
 */
class AngstromWrapper {
public:
  Scine::Utils::PositionCollection positions;

  AngstromWrapper() = default;
  explicit AngstromWrapper(unsigned N);
  explicit AngstromWrapper(
    Scine::Utils::PositionCollection pos,
    LengthUnit lengthUnit = LengthUnit::Bohr
  );

  /*!
   * @brief Fetch a bohr representation of the wrapped positions
   *
   * @warning After calling this function, you should not reuse the
   * corresponding instance, as the underlying positions have been converted to
   * bohr.
   */
  Scine::Utils::PositionCollection getBohr();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
  bool _invalidated = false;
};

} // namespace molassmbler

} // namespace Scine

#endif
