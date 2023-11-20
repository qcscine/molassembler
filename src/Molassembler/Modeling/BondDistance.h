/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Bond distance modelling functions.
 */

#ifndef INCLUDE_MOLASSEMBLER_BOND_DISTANCE_H
#define INCLUDE_MOLASSEMBLER_BOND_DISTANCE_H

#include "Molassembler/Modeling/AtomInfo.h"
#include "Molassembler/Types.h"

namespace Scine {
namespace Molassembler {
namespace Bond {

//! Bond order definition for bond types as defined in common_typedefs
static constexpr std::array<double, 7> bondOrderMap {{
  1, 2, 3, 4, 5, 6, 0.5
}};

//! UFF bond distance correction constant lambda
constexpr double bondOrderCorrectionLambda = 0.1332;

/*! @brief Calculates bond distance as modeled by UFF
 *
 * @complexity{@math{\Theta(1)}}
 */
double calculateBondDistance(
  Utils::ElementType a,
  Utils::ElementType b,
  BondType bondType
);

/*! @brief Calculates bond order as modeled by UFF
 *
 * @complexity{@math{\Theta(1)}}
 */
double calculateBondOrder(
  Utils::ElementType a,
  Utils::ElementType b,
  double distance
);

} // namespace Bond
} // namespace Molassembler
} // namespace Scine

#endif
