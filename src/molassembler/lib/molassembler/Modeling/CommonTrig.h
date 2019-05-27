/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Trigonometric stubs
 *
 * Contains some common trigonometric functionality. Does not contain custom
 * implementations of the basic trigonometric functions, but rather common
 * combinations of them in specific contexts.
 */

#ifndef INCLUDE_MOLASSEMBLER_COMMON_TRIG_H
#define INCLUDE_MOLASSEMBLER_COMMON_TRIG_H

#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/Modeling/BondDistance.h"

#include <cmath>

namespace Scine {

namespace molassembler {

namespace CommonTrig {

using ValueBounds = DistanceGeometry::ValueBounds;

/*!
 * Calculates the law of cosines, returning c:
 * c² = a² + b² - 2ab cos φ
 */
template<typename T>
T lawOfCosines(
  const T a,
  const T b,
  const T phiRadians
) {
  return std::sqrt(
    a * a
    + b * b
    - 2 * a * b * std::cos(phiRadians)
  );
}

/*!
 * Calculates the angle between a and b in a triangle of three side lengths a,
 * b and c.
 */
template<typename T>
T lawOfCosinesAngle(
  const T a,
  const T b,
  const T c
) {
  double ratio = (
    a * a
    + b * b
    - c * c
  ) / (
    2 * a * b
  );

  /* Accept ratios a little above one for the sake of numerical inaccuracies
   * when the triangle almost doesn't exist
   */
  if(ratio > 1 && std::fabs(ratio - 1) <= 1e-10) {
    return 0.0;
  }

  return std::acos(ratio);
}

/*!
 * Calculates the law of sines, returning the angle β.
 * From: (sin α) / a = (sin β) / b
 * -> β = arcsin( (b sin α) / a )
 */
template<typename T>
T lawOfSinesAngle(
  const T a,
  const T alphaRadians,
  const T b
) {
  return std::asin(
    b * std::sin(alphaRadians) / a
  );
}

/*!
 * Calculates the 1-4 length in a dihedral, which is defined by the three
 * distances a, b, c, the two angles at ab and at bc and the dihedral around the
 * intermediate bond b.
 */
double dihedralLength(
  double a,
  double b,
  double c,
  double alpha,
  double beta,
  double dihedral
);

/**
 * @brief Calculates the bounds on the dihedral length given bounds for all
 *   constituting variables
 *
 * @param aBounds The bond length bounds between atoms i and j
 * @param bBounds The bond length bounds between atoms j and k
 * @param cBounds The bond length bounds between atoms k and l
 * @param alphaBounds The angle bounds for the sequence i-j-k
 * @param betaBounds The angle bounds for the sequence j-k-l
 * @param dihedralBounds The dihedral bounds for the sequence i-j-k-l
 *
 * @return Bounds on the dihedral length
 */
ValueBounds dihedralLengthBounds(
  const ValueBounds& aBounds,
  const ValueBounds& bBounds,
  const ValueBounds& cBounds,
  const ValueBounds& alphaBounds,
  const ValueBounds& betaBounds,
  const ValueBounds& dihedralBounds
);

} // namespace CommonTrig

} // namespace molassembler

} // namespace Scine

#endif
