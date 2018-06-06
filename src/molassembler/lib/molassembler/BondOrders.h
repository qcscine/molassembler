#ifndef INCLUDE_MOLASSEMBLER_BOND_ORDERS_H
#define INCLUDE_MOLASSEMBLER_BOND_ORDERS_H

#include "Delib/BondOrderCollection.h"
#include "Delib/ElementTypeCollection.h"

#include "detail/AngstromWrapper.h"

namespace molassembler {

//! Calculates a bond order collection via UFF-like bond distance modelling
Delib::BondOrderCollection uffBondOrders(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper
);

} // namespace molassembler

#endif
