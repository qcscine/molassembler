// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_BOND_ORDERS_H
#define INCLUDE_MOLASSEMBLER_BOND_ORDERS_H

#include "Utils/ElementTypes.h"
#include <vector>

// Forward-declarations
namespace Scine {
namespace Utils {
class BondOrderCollection;
using ElementTypeCollection = std::vector<ElementType>;
} // namespace Utils
} // namespace Scine

/*! @file
 *
 * @brief Interpret bond orders from 3D coordinates only
 */

namespace molassembler {

// Forward-declarations
class AngstromWrapper;

/*!
 * @brief Calculates a bond order collection via UFF-like bond distance modelling
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

#endif
