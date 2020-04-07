/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
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

/*! @brief Calculates bond distance as modelled by UFF
 *
 * @complexity{@math{\Theta(1)}}
 */
double calculateBondDistance(
  Utils::ElementType a,
  Utils::ElementType b,
  BondType bondType
);

/*! @brief Calculates bond distances as modelled by UFF
 *
 * @complexity{@math{\Theta(1)}}
 */
double calculateBondOrder(
  Utils::ElementType a,
  Utils::ElementType b,
  double distance
);

} // namespace Bond

} // namespace molassembler

} // namespace Scine

#endif
