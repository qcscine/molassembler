#include "CommonTrig.h"

namespace MoleculeManip {

namespace CommonTrig {


/* TODO maybe this belongs into a StereocentersHelpers.h
 *  along with many other such cases...
 */
double getRingOneFourDistance(
  const double& a,
  const double& b,
  const double& c,
  const double& insideAngleRadians
) {
  /*
   *            d
   *   1 - - - - - - - - 4
   *    \        _ .·°  /
   *   a \   _.·°      / c
   *      \α_________α/
   *      2     b     3
   *
   * 1-2: a
   * 2-3: b
   * 3-4: c
   * 2-4: x <- law of cosines (b, c, α)
   * angle β (3-2-4) <- law of sines (α, x, c)
   * 1-4: d <- law of cosines (a, x, α-β)
   *
   * Although this figure suggests it, b is not parallel to d!
   */
  const double x = CommonTrig::lawOfCosines(
    b,
    c,
    insideAngleRadians
  );
  return CommonTrig::lawOfCosines(
    a,
    x,
    insideAngleRadians - CommonTrig::lawOfSinesAngle(
      x,
      insideAngleRadians,
      c
    )
  );
}

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
    + 4 * b * b
    + c * c
    - 4 * a * b * cos(abAngle)
    + 2 * a * c * (
      sin(abAngle) * sin(bcAngle) * cos(dihedral) 
      - cos(abAngle) * cos(bcAngle)
    )
    + 4 * b * c * cos(bcAngle)
  );
}

} // eo CommonTrig

}
