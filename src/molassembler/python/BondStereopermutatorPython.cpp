/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "OptionalPython.h"
#include "pybind11/operators.h"
#include "pybind11/stl.h"

#include "molassembler/AtomStereopermutator.h"
#include "molassembler/BondStereopermutator.h"

void init_bond_stereopermutator(pybind11::module& m) {
  using namespace Scine::molassembler;

  pybind11::class_<BondStereopermutator> bondStereopermutator(
    m,
    "BondStereopermutator",
    R"delim(
      Handles specific relative arrangements of two atom stereopermutators
      joined by a bond. This includes, importantly, E/Z stereocenters at double
      bonds.

      >>> # The bond stereopermutator in but-2z-ene
      >>> import molassembler as masm
      >>> z_butene = masm.io.experimental.from_smiles("C/C=C\C")
      >>> bond_index = masm.BondIndex(1, 2)
      >>> assert z_butene.graph.bond_type(bond_index) = masm.BondType.Double
      >>> permutator = masm.stereopermutators.option(bond_index)
      >>> assert permutator is not None
      >>> permutator.assigned is not None
      True
      >>> permutator.num_assignments
      2
    )delim"
  );

  pybind11::enum_<BondStereopermutator::FittingMode> fittingMode(
    bondStereopermutator,
    "FittingMode",
    "Differentiates how viable assignments are chosen during fitting"
  );
  fittingMode.value(
    "Thresholded",
    BondStereopermutator::FittingMode::Thresholded,
    "Positions must be close to the idealized assignment geometry"
  );
  fittingMode.value(
    "Nearest",
    BondStereopermutator::FittingMode::Nearest,
    "The assignment closest to the idealized geometry is chosen"
  );


  bondStereopermutator.def(
    "assigned",
    &BondStereopermutator::assigned,
    R"delim(
      Returns an integer indicating the assignment of the stereopermutator or
      ``None`` if the stereopermutator is unassigned.

      >>> # An unassigned bond stereopermutator
      >>> import molassembler as masm
      >>> butene = masm.io.experimental.from_smiles("CC=CC")
      >>> bond_index = masm.BondIndex(1, 2)
      >>> assert z_butene.graph.bond_type(bond_index) = masm.BondType.Double
      >>> permutator = masm.stereopermutators.option(bond_index)
      >>> assert permutator is not None
      >>> permutator.assigned is None
      True
      >>> permutator.num_assignments
      2
    )delim"
  );

  bondStereopermutator.def(
    "index_of_permutation",
    &BondStereopermutator::indexOfPermutation,
    R"delim(
      Returns an integer indicating the index of permutation if the
      stereopermutator is assigned or ``None`` if the stereopermutator is
      unassigned.

      >>> # A case in which the number of abstract and feasible permutations
      >>> # differ: bond stereopermutators in small cycles (<= 6)
      >>> import molassembler as masm
      >>> benzene = masm.io.experimental.from_smiles("C1=CC=CC=C1")
      >>> permutators = benzene.stereopermutators.bond_stereopermutators()
      >>> has_two_stereopermutations = lambda p: p.num_stereopermutations == 2
      >>> has_one_assignment = lamda p: p.num_assignments == 1
      >>> all(map(has_two_stereopermutations, permutators))
      True
      >>> all(map(has_one_assignment, permutators))
      True
    )delim"
  );

  bondStereopermutator.def(
    "num_assignments",
    &BondStereopermutator::numAssignments,
    R"delim(
      Returns the number of assignments. Valid assignment indices range from 0
      to this number minus one.
    )delim"
  );

  bondStereopermutator.def(
    "num_stereopermutations",
    &BondStereopermutator::numStereopermutations,
    "Returns the number of stereopermutations."
  );

  bondStereopermutator.def(
    "edge",
    &BondStereopermutator::edge,
    "Returns the edge this stereopermutator is placed on."
  );

  bondStereopermutator.def(
    "dihedral",
    &BondStereopermutator::dihedral,
    pybind11::arg("stereopermutator_a"),
    pybind11::arg("site_index_a"),
    pybind11::arg("stereopermutator_b"),
    pybind11::arg("site_index_b"),
    R"delim(
      Returns the dihedral angle between two sites of the constituting atom
      stereopermutators in radians

      You can glean site indices from the individual constituting atom
      stereopermutators' rankings.

      :param stereopermutator_a: One constituting :class:`AtomStereopermutator`
      :param site_index_a: The site index of ``stereopermutator_a`` starting the dihedral sequence
      :param stereopermutator_b: The other constituting :class:`AtomStereopermutator`
      :param site_index_b: The site index of ``stereopermutator_b`` ending the dihedral sequence
    )delim"
  );

  bondStereopermutator.def(pybind11::self == pybind11::self);
  bondStereopermutator.def(pybind11::self != pybind11::self);

  bondStereopermutator.def(
    "__repr__",
    [](const BondStereopermutator& permutator) {
      return permutator.info();
    }
  );
}
