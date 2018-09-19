// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Modeling/CommonTrig.h"

namespace molassembler {

namespace CommonTrig {

double dihedralLength(
  const double a,
  const double b,
  const double c,
  const double abAngle,
  const double bcAngle,
  const double dihedral
) {
  return sqrt(
    a * a
    + b * b
    + c * c
    + 2 * (
      - a * b * std::cos(abAngle)
      + a * c * (
        std::cos(abAngle) * std::cos(bcAngle)
        - std::sin(abAngle) * std::sin(bcAngle) * std::cos(dihedral)
      )
      - b * c * std::cos(bcAngle)
    )
  );
}

} // namespace CommonTrig

} // namespace molassembler
