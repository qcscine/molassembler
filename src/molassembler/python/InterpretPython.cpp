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
    "Mapping of atom indices from the original positional information to which "
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
    R"delim(
      Interpret molecules from element types, positional information and bond orders

      Attempts to interpret (possibly multiple) Molecules from element types,
      positional information and a bond order collection. Bond orders are
      discretized into bond types. Connected components within the space are
      identified and individually instantiated into Molecules. The
      instantiation of BondStereopermutators in the Molecules can be limited to
      edges whose bond order exceeds a particular value.

      :param atom_collection: Element types and positional information in Bohr units
      :param bond_orders: Fractional bond orders
      :param discretization: How bond fractional orders are to be discretized
      :param stereopermutator_bond_order_threshold: If specified, limits the
        instantiation of BondStereopermutators onto edges whose fractional bond
        orders exceed the provided threshold. If None, BondStereopermutators
        are instantiated at all bonds.
      :raises ValueError: If the number of particles in the atom collection and
        bond order collections do not match
    )delim"
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
    R"delim(
      Interpret molecules from element types and positional information. Bond
      orders are calculated with UFF parameters.

      Attempts to interpret (possibly multiple) Molecules from element types
      and positional information. Bond orders are calculated from atom-pairwise
      spatial distances using UFF parameters. The bond orders are then
      discretized into bond types. Connected components within the space are
      identified and individually instantiated into Molecules. The
      instantiation behavior of BondStereopermutators in the Molecules can be
      limited to edges whose bond order exceeds a particular value.

      :param atom_collection: Element types and positional information in Bohr units
      :param discretization: How bond fractional orders are to be discretized
      :param stereopermutator_bond_order_threshold: If specified, limits the
        instantiation of BondStereopermutators onto edges whose fractional bond orders
        exceed the provided threshold. If None, BondStereopermutators are
        instantiated at all bonds.
    )delim"
  );

  m.def(
    "apply_interpretation_map",
    &applyInterpretationMap,
    pybind11::arg("interpret_result"),
    pybind11::arg("atom_collection"),
    R"delim(
      Splits an atom collection just like an interpret split the positions
      into multiple molecules

      :param interpret_result: The result of an interpret call on the same atom collection
      :param atom_collection: The atom collection to split
      :rtype: List of atom collections
    )delim"
  );
}
