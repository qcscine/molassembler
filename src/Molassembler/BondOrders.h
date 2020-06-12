/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Interpret bond orders from 3D coordinates only
 */

#ifndef INCLUDE_MOLASSEMBLER_BOND_ORDERS_H
#define INCLUDE_MOLASSEMBLER_BOND_ORDERS_H

#include "Molassembler/Export.h"
#include "Utils/Geometry/ElementTypes.h"
#include <vector>

namespace Scine {
namespace Utils {
class BondOrderCollection;
using ElementTypeCollection = std::vector<ElementType>;
} // namespace Utils

namespace Molassembler {

// Forward-declarations
class AngstromPositions;

/*! @brief Calculates a floating-point bond order collection via UFF-like bond distance modelling
 *
 * @complexity{@math{\Theta(N^2)}}
 * @throws std::logic_error If interpreted fractional bond orders are greater
 *   than 6.5.  In these cases, the structure is most likely unreasonable.
 * @warning UFF parameter bond order calculation is a very primitive
 *   approximation and carries a high risk of misinterpretation
 */
MASM_EXPORT Utils::BondOrderCollection uffBondOrders(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromPositions
);

/*! @brief Calculates a binary (single or none) bond order collection via covalent radii
 *
 * @complexity{@math{\Theta(N^2)}}
 * @note This is just a convenience forwarder for
 * Utils::BondDetector::detectBonds for different types
 */
MASM_EXPORT Utils::BondOrderCollection covalentRadiiBondOrders(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromPositions
);

} // namespace Molassembler
} // namespace Scine

#endif
