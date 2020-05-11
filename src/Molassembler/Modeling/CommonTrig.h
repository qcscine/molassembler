/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Trigonometric stubs
 *
 * Contains some common trigonometric functionality. Does not contain custom
 * implementations of the basic trigonometric functions, but rather common
 * combinations of them in specific contexts.
 */

#ifndef INCLUDE_MOLASSEMBLER_COMMON_TRIG_H
#define INCLUDE_MOLASSEMBLER_COMMON_TRIG_H

#include "Molassembler/DistanceGeometry/ValueBounds.h"
#include "Molassembler/Modeling/BondDistance.h"

#include <cmath>

namespace Scine {

namespace Molassembler {

namespace CommonTrig {

using ValueBounds = DistanceGeometry::ValueBounds;

/*! @brief Calculates the law of cosines
 *
 * Calculates the law of cosines, returning @math{c} from
 * @math{c^2 = a^2 + b^2 - 2ab \cos \varphi}
 *
 * @complexity{@math{\Theta(1)}}
 */
template<typename T>
T lawOfCosines(const T a, const T b, const T phiRadians) {
  return std::sqrt(
    a * a
    + b * b
    - 2 * a * b * std::cos(phiRadians)
  );
}

/*! @brief Calculates angle between a and b in a triangle of three side lengths
 *
 * @complexity{@math{\Theta(1)}}
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

/*! @brief Calculate the angle beta using the law of sines
 *
 * Calculates the law of sines, returning the angle @math{\beta}.
 * From: @math{\frac{\sin \alpha}{a} = \frac{\sin \beta}{b}}
 * To: @math{\beta = \mathrm{arcsin}(b \frac{\sin \alpha}{a})}
 *
 * @complexity{@math{\Theta(1)}}
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

/*! @brief Calculates the distance between the 1st and 4th position in a
 *   dihedral
 *
 * Calculates the 1-4 length in a dihedral, which is defined by the three
 * distances a, b, c, the two angles at ab and at bc and the dihedral around the
 * intermediate bond b.
 *
 * @complexity{@math{\Theta(1)}}
 */
double dihedralLength(
  double a,
  double b,
  double c,
  double alpha,
  double beta,
  double dihedral
);

/** @brief Calculates the bounds on the dihedral length given bounds for all
 *   constituting variables
 *
 * @complexity{Unclear as this is minimized, not calculated directly.}
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

} // namespace Molassembler

} // namespace Scine

#endif
