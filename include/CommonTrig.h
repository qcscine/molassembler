#ifndef INCLUDE_COMMON_TRIG_H
#define INCLUDE_COMMON_TRIG_H

#include <cmath>
#include "ElementTypes.h" // Delib

#include "BondDistance.h"

namespace MoleculeManip {

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

/* TODO maybe this belongs into a StereocentersHelpers.h
 *  along with many other such cases...
 */
double getRingOneFourDistance(
  const double& a,
  const double& b,
  const double& c,
  const double& insideAngleRadians
);

} // eo CommonTrig

}

#endif
