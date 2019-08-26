/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Interpret bond orders from 3D coordinates only
 */

#ifndef INCLUDE_MOLASSEMBLER_BOND_ORDERS_H
#define INCLUDE_MOLASSEMBLER_BOND_ORDERS_H

#include "Utils/Geometry/ElementTypes.h"
#include <vector>

// Forward-declarations
namespace Scine {
namespace Utils {
class BondOrderCollection;
using ElementTypeCollection = std::vector<ElementType>;
} // namespace Utils
} // namespace Scine

namespace Scine {

namespace molassembler {

// Forward-declarations
class AngstromWrapper;

/*! @brief Calculates a bond order collection via UFF-like bond distance modelling
 *
 * @complexity{@math{\Theta(N^2)}}
 * @throws std::logic_error If interpreted fractional bond orders are greater
 *   than 6.5.  In these cases, the structure is most likely unreasonable.
 * @warning UFF parameter bond order calculation is a very primitive
 *   approximation and carries a high risk of misinterpretation
 */
Scine::Utils::BondOrderCollection uffBondOrders(
  const Scine::Utils::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper
);

} // namespace molassembler

} // namespace Scine

#endif
