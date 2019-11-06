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
    "Handles specific relative arrangements of two atom stereopermutators "
    "joined by a bond"
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
    "Returns an integer indicating the assignment of the stereopermutator or "
    "None if the stereopermutator is unassigned."
  );

  bondStereopermutator.def(
    "index_of_permutation",
    &BondStereopermutator::indexOfPermutation,
    "Returns an integer indicating the index of permutation if the "
    "stereopermutator is assigned or None if the stereopermutator is unassigned."
  );

  bondStereopermutator.def(
    "num_assignments",
    &BondStereopermutator::numAssignments,
    "Returns the number of assignments. Valid assignment indices range from 0 "
    "to this number minus one."
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
      Returns the dihedral angle between two sites of the constituting atom stereopermutators

      You can glean site indices from the individual constituting atom
      stereopermutators' rankings.

      :param stereopermutator_a: One constituting atom stereopermutator
      :param site_index_a: The site index of stereopermutator_a starting the dihedral sequence
      :param stereopermutator_b: The other constituting atom stereopermutator
      :param site_index_b: The site index of stereopermutator_b ending the dihedral sequence
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
