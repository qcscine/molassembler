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
    "Returns an optional type indicating whether the stereopermutator is "
    "assigned, and if so, which index of assignment."
  );

  bondStereopermutator.def(
    "index_of_permutation",
    &BondStereopermutator::indexOfPermutation,
    "Returns an optional type indicating whether the stereopermutator is "
    "assigned, and if so, which index of permutation."
  );

  bondStereopermutator.def(
    "num_assignments",
    &BondStereopermutator::numAssignments,
    "Returns the number of possible assignments"
  );

  bondStereopermutator.def(
    "num_stereopermutations",
    &BondStereopermutator::numStereopermutations,
    "Returns the number of stereopermutations"
  );

  bondStereopermutator.def(
    "edge",
    &BondStereopermutator::edge,
    "Returns the edge this stereopermutator is placed on"
  );

  bondStereopermutator.def(
    "dihedral",
    &BondStereopermutator::dihedral,
    pybind11::arg("stereopermutator_a"),
    pybind11::arg("site_index_a"),
    pybind11::arg("stereopermutator_b"),
    pybind11::arg("site_index_b"),
    "Returns the dihedral angle between two sites of the constituting atom stereopermutators"
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
