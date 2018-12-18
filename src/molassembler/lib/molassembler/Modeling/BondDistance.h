/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Bond distance modelling functions.
 */

#ifndef INCLUDE_MOLASSEMBLER_BOND_DISTANCE_H
#define INCLUDE_MOLASSEMBLER_BOND_DISTANCE_H

#include "molassembler/Modeling/AtomInfo.h"
#include "molassembler/Types.h"

namespace Scine {

namespace molassembler {

namespace Bond {

//! Bond order definition for bond types as defined in common_typedefs
static constexpr std::array<double, 7> bondOrderMap {{
  1, 2, 3, 4, 5, 6, 0.5
}};

//! UFF bond distance correction constant lambda
constexpr double bondOrderCorrectionLambda = 0.1332;

//! Calculates bond distance as modelled by UFF
double calculateBondDistance(
  Delib::ElementType a,
  Delib::ElementType b,
  BondType bondType
);

//! Calculates bond distances as modelled by UFF
double calculateBondOrder(
  Delib::ElementType a,
  Delib::ElementType b,
  double distance
);

} // namespace Bond

} // namespace molassembler

} // namespace Scine

#endif
