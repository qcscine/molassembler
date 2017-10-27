#include "CommonTrig.h"

namespace MoleculeManip {

namespace CommonTrig {

double dihedralLength(
  const double& a,
  const double& b,
  const double& c,
  const double& abAngle,
  const double& bcAngle,
  const double& dihedral
) {
  return sqrt(
    a * a
    + b * b
    + c * c
    + 2 * (
      - a * b * cos(abAngle)
      + a * c * (
        cos(abAngle) * cos(bcAngle)
        - sin(abAngle) * sin(bcAngle) * cos(dihedral)
      )
      - b * c * cos(bcAngle)
    )
  );
}

} // namespace CommonTrig

} // namespace MoleculeManip
