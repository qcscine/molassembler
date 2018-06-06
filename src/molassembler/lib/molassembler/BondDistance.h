#ifndef INCLUDE_BOND_DISTANCE_H
#define INCLUDE_BOND_DISTANCE_H

#include "AtomInfo.h"

#include "detail/SharedTypes.h"

#include <cmath>
#include <array>

/*! @file
 *
 * Bond distance modelling functions.
 */

namespace molassembler {

namespace Bond {

//! Bond order definition for bond types as defined in common_typedefs
static constexpr std::array<double, 8> bondOrderMap {{
  1, 2, 3, 4, 5, 6, 1.5, 0.5
}};

//! UFF bond distance correction constant lambda
constexpr double bondOrderCorrectionLambda = 0.1332;

//! Calculates bond distance as modelled by UFF
double calculateBondDistance(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
);

//! Calculates bond distances as modelled by UFF
double calculateBondOrder(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const double& distance
);

} // namespace Bond

} // namespace molassembler

#endif
