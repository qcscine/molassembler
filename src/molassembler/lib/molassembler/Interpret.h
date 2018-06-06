#ifndef INCLUDE_MOLASSEMBLER_INTERPRET_H
#define INCLUDE_MOLASSEMBLER_INTERPRET_H

#include "Molecule.h"

#include "Delib/AtomCollection.h"
#include "Delib/BondOrderCollection.h"

namespace molassembler {

//! Specify the algorithm used to discretize floating-point bond orders into bond types
enum class BondDiscretizationOption : unsigned {
  //! All bond orders >= 0.5 are considered single bonds
  Binary,
  //! Bond orders are rounded to the nearest integer
  UFF
};

//! Result type of an interpret call.
struct InterpretResult {
  //! The list of individual molecules found in the 3D information
  std::vector<Molecule> molecules;
  //! A map from an index within the AtomCollection to its index in molecules member
  std::vector<unsigned> componentMap;
};

InterpretResult interpret(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const Delib::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization = BondDiscretizationOption::Binary
);

InterpretResult interpret(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const BondDiscretizationOption discretization = BondDiscretizationOption::Binary
);

/*! Interpret molecules in 3D information and a bond order collection.
 *
 * The graph is inferred via bond discretization from the bond order
 * collection.
 *
 * @note Assumes that the provided atom collection's positions are in
 * Bohr units.
 */
InterpretResult interpret(
  const Delib::AtomCollection& atomCollection,
  const Delib::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization = BondDiscretizationOption::Binary
);

/*! Interpret molecules in 3D information.
 *
 * The graph is inferred via bond discretization from pairwise atom distances.
 *
 * @note Assumes that the provided atom collection's positions are in
 * Bohr units.
 */
InterpretResult interpret(
  const Delib::AtomCollection& atomCollection,
  const BondDiscretizationOption discretization = BondDiscretizationOption::Binary
);

} // namespace molassembler

#endif
