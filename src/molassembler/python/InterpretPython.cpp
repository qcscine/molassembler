/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "OptionalPython.h"

#include "molassembler/Interpret.h"
#include "molassembler/Molecule.h"

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Geometry/AtomCollection.h"

void init_interpret(pybind11::module& m) {
  using namespace Scine::molassembler;
  using namespace Scine::Utils;

  pybind11::enum_<BondDiscretizationOption>(
    m,
    "BondDiscretization",
    "Specifies the algorithm used to discretize floating-point bond orders into "
    "discrete bond types."
  ).value("Binary", BondDiscretizationOption::Binary, "All bond orders >= 0.5 are considered single bonds")
    .value("RoundToNearest", BondDiscretizationOption::RoundToNearest, "Round bond orders to nearest integer");

  pybind11::class_<InterpretResult> interpretResult(
    m,
    "InterpretResult",
    "Result type of an interpret call."
  );

  interpretResult.def_readwrite(
    "molecules",
    &InterpretResult::molecules,
    "Individual molecules found in the 3D information"
  );

  interpretResult.def_readwrite(
    "component_map",
    &InterpretResult::componentMap,
    "Mapping of indices from the original positional information to which "
    "molecule it is now part"
  );

  m.def(
    "interpret",
    pybind11::overload_cast<
      const AtomCollection&,
      const BondOrderCollection&,
      BondDiscretizationOption,
      const boost::optional<double>&
    >(&interpret),
    pybind11::arg("atom_collection"),
    pybind11::arg("bond_orders"),
    pybind11::arg("discretization"),
    pybind11::arg("stereopermutator_bond_order_threshold") = 1.4,
    "Interpret molecules from element types, positional information and bond orders"
  );

  m.def(
    "interpret",
    pybind11::overload_cast<
      const AtomCollection&,
      BondDiscretizationOption,
      const boost::optional<double>&
    >(&interpret),
    pybind11::arg("atom_collection"),
    pybind11::arg("discretization"),
    pybind11::arg("stereopermutator_bond_order_threshold") = 1.4,
    "Interpret molecules from element types and positional information. Bond "
    "orders are then interpreted using UFF parameters."
  );
}
