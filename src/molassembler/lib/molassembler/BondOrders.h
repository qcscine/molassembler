/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Interpret bond orders from 3D coordinates only
 */

#ifndef INCLUDE_MOLASSEMBLER_BOND_ORDERS_H
#define INCLUDE_MOLASSEMBLER_BOND_ORDERS_H

// Forward-declarations
namespace Delib {
class BondOrderCollection;
class ElementTypeCollection;
} // namespace Delib

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
Delib::BondOrderCollection uffBondOrders(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper
);

} // namespace molassembler

#endif
