/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

#include "Molassembler/Interpret.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Graph.h"

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Geometry/AtomCollection.h"

void init_interpret(pybind11::module& m) {
  using namespace Scine::Molassembler;
  using namespace Scine::Utils;

  auto interpretSubmodule = m.def_submodule("interpret");
  interpretSubmodule.doc() = R"delim(
    Submodule with freestanding functions yielding :class:`~scine_molassembler.Molecule` or
    :class:`~scine_molassembler.Graph` instances from Cartesian coordinate data (and optionally
    bond order data).

    **Bond discretization**

    In the discretization of fractional bond orders to classic integer internal
    bond types (e.g. single, double, etc.), there are two options. You can
    choose to round bond orders to the nearest integer, but this is
    particularly error prone for particularly weakly-bound metal ligands
    (around 0.5) and aromatic bonds (around 1.5). For instance, two adjacent
    aromatic bonds that both show a fractional bond order around 1.5 may be
    randomly rounded up or down depending on the bond order generation method
    or its particular conformation. This can cause unexpected ranking
    inequivalency / equivalency artifacts. If you expect there to be conjugated
    systems or transition metals in your set of interpreted molecules,
    discretizing bond orders in this fashion is currently disadvised.

    It can instead be preferable to discretize bond orders in a purely binary
    manner, i.e. bond orders are interpreted as a single bond if the fractional
    bond order is is more than or equal to 0.5. Double bond stereocenters (i.e.
    in organic molecules E/Z stereocenters) are still interpreted from
    coordinate information despite the main bond type discretized to a single
    bond.
  )delim";

  pybind11::enum_<Interpret::BondDiscretizationOption>(
    interpretSubmodule,
    "BondDiscretization",
    R"delim(
      Specifies the algorithm used to discretize floating-point bond orders into
      discrete bond types.
    )delim"
  ).value("Binary", Interpret::BondDiscretizationOption::Binary, "All bond orders >= 0.5 are considered single bonds")
    .value("RoundToNearest", Interpret::BondDiscretizationOption::RoundToNearest, "Round bond orders to nearest integer");

  pybind11::class_<Interpret::MoleculesResult> interpretResult(
    interpretSubmodule,
    "MoleculesResult",
    "Result type of a molceule interpret call."
  );

  interpretResult.def_readwrite(
    "molecules",
    &Interpret::MoleculesResult::molecules,
    R"delim(
      Individual molecules found in the 3D information.

      :rtype: ``List`` of :class:`~scine_molassembler.Molecule`
    )delim"
  );

  interpretResult.def_readwrite(
    "component_map",
    &Interpret::MoleculesResult::componentMap,
    R"delim(
      Mapping of atom indices from the original positional information to which
      molecule it is now part.

      :rtype: ``List[int]``
    )delim"
  );

  interpretSubmodule.def(
    "molecules",
    pybind11::overload_cast<
      const AtomCollection&,
      const BondOrderCollection&,
      Interpret::BondDiscretizationOption,
      const boost::optional<double>&
    >(&Interpret::molecules),
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
        orders exceed the provided threshold. If ``None``, BondStereopermutators
        are instantiated at all bonds.
      :raises ValueError: If the number of particles in the atom collection and
        bond order collections do not match

      >>> import scine_utilities as utils
      >>> import numpy as np
      >>> elements = [utils.ElementType.H] * 4
      >>> positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.71, 0.0], [2.0, 2.0, 2.0], [2.0, 2.71, 2.0]])
      >>> atoms = utils.AtomCollection(elements, positions)
      >>> bond_orders = utils.BondOrderCollection(4)
      >>> bond_orders.set_order(0, 1, 1.0)
      >>> bond_orders.set_order(2, 3, 1.0)
      >>> discretization = BondDiscretization.RoundToNearest
      >>> result = interpret(atoms, bond_orders, discretization)
      >>> assert len(result.molecules) == 2
      >>> hydrogen = Molecule()
      >>> assert all([m == hydrogen for m in result.molecules])
      >>> assert result.component_map == [0, 0, 1, 1]
    )delim"
  );

  interpretSubmodule.def(
    "interpret",
    pybind11::overload_cast<
      const AtomCollection&,
      Interpret::BondDiscretizationOption,
      const boost::optional<double>&
    >(&Interpret::molecules),
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
        exceed the provided threshold. If ``None``, BondStereopermutators are
        instantiated at all bonds.
    )delim"
  );

  interpretSubmodule.def(
    "apply_interpretation_map",
    &Interpret::applyInterpretationMap,
    pybind11::arg("component_map"),
    pybind11::arg("atom_collection"),
    R"delim(
      Splits an atom collection just like an interpret split the positions
      into multiple molecules

      :param component_map: The component map of an interpret call on the same atom collection
      :param atom_collection: The atom collection to split
      :rtype: List of atom collections
    )delim"
  );

  pybind11::class_<Interpret::GraphsResult> graphsResult(
    interpretSubmodule,
    "GraphsResult",
    "Result type of a graph interpret call."
  );

  graphsResult.def_readwrite(
    "graphs",
    &Interpret::GraphsResult::graphs,
    R"delim(
      Individual graphs found in the 3D information.

      :rtype: ``List`` of :class:`~scine_molassembler.Graph`
    )delim"
  );

  graphsResult.def_readwrite(
    "component_map",
    &Interpret::GraphsResult::componentMap,
    R"delim(
      Mapping of atom indices from the original positional information to which
      molecule it is now part.

      :rtype: ``List[int]``
    )delim"
  );

  interpretSubmodule.def(
    "graphs",
    pybind11::overload_cast<
      const AtomCollection&,
      const BondOrderCollection&,
      Interpret::BondDiscretizationOption
    >(&Interpret::graphs),
    pybind11::arg("atom_collection"),
    pybind11::arg("bond_orders"),
    pybind11::arg("discretization"),
    R"delim(
      Interpret graphs from element types, positional information and bond orders

      Attempts to interpret (possibly multiple) graphs from element types,
      positional information and a bond order collection. Bond orders are
      discretized into bond types. Connected components within the space are
      identified and individually instantiated into graphs.

      :param atom_collection: Element types and positional information in Bohr units
      :param bond_orders: Fractional bond orders
      :param discretization: How bond fractional orders are to be discretized
      :raises ValueError: If the number of particles in the atom collection and
        bond order collections do not match
    )delim"
  );

  interpretSubmodule.def(
    "invert_component_map",
    &Interpret::invertComponentMap,
    pybind11::arg("component_map"),
    R"delim(
      Inverts a component mapping.

      Allows direct determination of the original index of a component's atom
      within the coordinate set used in interpretation.

      :param component_map: A component map as yielded in :class:`GraphsResult`
        or :class:`MoleculesResult`.
      :returns: A nested list that contains the original indices for each
        component.

      >>> component_map = [0, 1, 1, 0, 1] # 0->0, 1->1, 2->1, etc.
      >>> inv = invert_component_map(component_map)
      >>> inv
      [[0, 3], [1, 2, 4]]
      >>> def lookup(component, index):
      ...     return inv[component][index]
      >>> lookup(component=1, index=2)
      4
      >>> lookup(component=0, index=1)
      3
    )delim"
  );

  pybind11::class_<Interpret::FalsePositive> falsePositive(
    interpretSubmodule,
    "FalsePositive",
    "Tuple-imitating data struct for false positives in bond discretization"
  );
  falsePositive.def_readwrite("i", &Interpret::FalsePositive::i);
  falsePositive.def_readwrite("j", &Interpret::FalsePositive::j);
  falsePositive.def_readwrite("probability", &Interpret::FalsePositive::probability, "Probability that a bond is a false positive by some arbitrary measure. Normalized between 0 and 1");
  falsePositive.def("__getitem__", [](const Interpret::FalsePositive& fp, const unsigned i) -> boost::variant<unsigned, double> {
    if(i == 0) {
      return fp.i;
    }

    if(i == 1) {
      return fp.j;
    }

    if(i == 2) {
      return fp.probability;
    }

    throw std::out_of_range("Only three elements in this tuple-like object");
  });

  falsePositive.def("__repr__", [](const Interpret::FalsePositive& fp) -> std::string {
    return "(" + std::to_string(fp.i) + ", "+ std::to_string(fp.j) + ", " + std::to_string(fp.probability) + ")";
  });

  interpretSubmodule.def(
    "uncertain_bonds",
    &Interpret::uncertainBonds,
    pybind11::arg("atom_collection"),
    pybind11::arg("bond_collection"),
    R"delim(
      Lists bonds with uncertain shape classifications at both ends

      :returns: a list of :class:`FalsePositive` objects

      .. warning::

         Do not alter both bonds if there is a bond pair that have overlapping
         indices. If suggested bonds overlap, remove only that bond with the
         higher probability
    )delim"
  );

  interpretSubmodule.def(
    "bad_haptic_ligand_bonds",
    &Interpret::badHapticLigandBonds,
    pybind11::arg("atom_collection"),
    pybind11::arg("bond_collection"),
    R"delim(
      Suggest false positive haptic ligand bonds

      Generates a plane of best fit for each haptic ligand in the interpreted
      graphs. If the angle of the normal of this plane to the axis defined by the
      central atom and the site centroid is more than 30 degrees, tries to name a
      single bond whose removal improves the interpretation.

      :returns: a list of :class:`FalsePositive` objects

      .. note::

         Suggested bonds can disconnect haptic sites. When making changes to a
         bond order matrix based on suggestions from this function, apply them
         one at a time based on the highest probability received. Additionally,
         if multiple bonds must be removed to make a haptic ligand
         geometrically reasonable, you will need to iteratively call this
         function and alter suggested bond orders.
    )delim"
  );

  interpretSubmodule.def(
    "remove_false_positives",
    &Interpret::removeFalsePositives,
    pybind11::arg("atoms"),
    pybind11::arg("bonds"),
    "Iteratively removes bonds reported by false positive detection functions"
  );
}
