#ifndef INCLUDE_COMMON_TRIG_H
#define INCLUDE_COMMON_TRIG_H

#include "BondDistance.h"

/*! @file
 *
 * Contains some common trigonometric functionality. Does not contain custom
 * implementations of the basic trigonometric functions, but rather common
 * combinations of them in specific contexts.
 *
 */

namespace molassembler {

namespace CommonTrig {

/*!
 * Calculates the law of cosines, returning c:
 * c² = a² + b² - 2ab cos φ
 */
template<typename T>
T lawOfCosines(
  const T& a,
  const T& b,
  const T& phiRadians
) {
  return sqrt(
    a * a
    + b * b
    - 2 * a * b * cos(phiRadians)
  );
}

/*!
 * Calculates the angle between a and b in a triangle of three side lengths a,
 * b and c.
 */
template<typename T>
T lawOfCosinesAngle(
  const T& a,
  const T& b,
  const T& c
) {
  return std::acos(
    (
      a * a
      + b * b
      - c * c
    ) / (
      2 * a * b
    )
  );
}

/*!
 * Calculates the law of sines, returning the angle β.
 * From: (sin α) / a = (sin β) / b
 * -> β = arcsin( (b sin α) / a )
 */
template<typename T>
T lawOfSinesAngle(
  const T& a,
  const T& alphaRadians,
  const T& b
) {
  return asin(
    b * sin(alphaRadians) / a
  );
}

/*!
 * Calculates the 1-4 length in a dihedral, which is defined by the three
 * distances a, b, c, the two angles at ab and at bc and the dihedral around the
 * intermediate bond b.
 */
double dihedralLength(
  const double& a,
  const double& b,
  const double& c,
  const double& abAngle,
  const double& bcAngle,
  const double& dihedral
);

} // namespace CommonTrig

} // namespace molassembler

#endif
