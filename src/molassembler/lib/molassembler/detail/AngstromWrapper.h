#ifndef INCLUDE_MOLASSEMBLER_ANGSTROM_POSITIONS_H
#define INCLUDE_MOLASSEMBLER_ANGSTROM_POSITIONS_H

#include "Delib/PositionCollection.h"

#include "molassembler/detail/SharedTypes.h"

/*!@file
 *
 * Provides a class that exists purely to strongly separate position collections
 * in bohr units and angstrom units in library interfaces
 */

namespace molassembler {

class AngstromWrapper {
public:
  Delib::PositionCollection positions;

  AngstromWrapper() = default;
  explicit AngstromWrapper(const unsigned N);
  explicit AngstromWrapper(
    Delib::PositionCollection pos,
    const LengthUnit lengthUnit = LengthUnit::Bohr
  );

  /*! Fetch a bohr representation of the wrapped positions
   *
   * \warning After calling this function, you should not reuse the
   * corresponding instance, as the underlying positions have been converted to
   * bohr.
   */
  Delib::PositionCollection getBohr();

private:
  bool _invalidated = false;
};

} // namespace molassmbler

#endif
