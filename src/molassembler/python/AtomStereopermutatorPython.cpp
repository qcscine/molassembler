/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "OptionalPython.h"
#include "pybind11/operators.h"
#include "pybind11/stl.h"

#include "molassembler/AtomStereopermutator.h"
#include "molassembler/RankingInformation.h"

void init_atom_stereopermutator(pybind11::module& m) {
  using namespace Scine::molassembler;
  pybind11::class_<AtomStereopermutator> atomStereopermutator(
    m,
    "AtomStereopermutator",
    "This class handles the permutation of ranked ligands around a central "
    "atom. It models its haptic ligands' binding sites and bridges in "
    "multidentate ligands in order to decide which stereopermutations are "
    "feasible. A stereopermutation may be infeasible, i.e. not realizable in "
    "three-dimensional space, if either haptic ligands would intersect due to "
    "too close ligand angles for their spatial heft, or if a multidentate "
    "ligand's bridge length between binding sites were too short to match the "
    "angle. The list of stereopermutations reduced by infeasible "
    "stereopermutations is then re-indexed and those indices into the list are "
    "called assignments.\n"
    "A Stereopermutator can be unassigned, i.e. the distinct stereopermutation "
    "that the substituents are can be indeterminate. If you choose to generate "
    "conformers for a molecule that includes unassigned stereopermutators, every "
    "conformer will choose an assignment from the pool of feasible assignments "
    "randomly, but consistent with relative statistical occurrence weights."
  );

  atomStereopermutator.def(
    "angle",
    &AtomStereopermutator::angle,
    pybind11::arg("ligand_index_i"),
    pybind11::arg("ligand_index_j"),
    "Fetches angle between substituent ligand indices in radians"
  );

  atomStereopermutator.def(
    "assigned",
    &AtomStereopermutator::assigned,
    "Returns an optional type indicating whether the stereopermutator is "
    "assigned, and if so, which index of assignment."
  );

  atomStereopermutator.def(
    "central_index",
    &AtomStereopermutator::centralIndex,
    "Returns the central atom this permutator is placed on"
  );

  atomStereopermutator.def(
    "index_of_permutation",
    &AtomStereopermutator::indexOfPermutation,
    "Returns an optional type indicating whether the stereopermutator is "
    "assigned, and if so, which index of permutation."
  );

  atomStereopermutator.def_property_readonly(
    "ranking",
    &AtomStereopermutator::getRanking,
    "Get the underlying ranking state of substituents"
  );

  atomStereopermutator.def_property_readonly(
    "symmetry",
    &AtomStereopermutator::getSymmetry,
    "Returns the underlying symmetry"
  );

  atomStereopermutator.def(
    "num_assignments",
    &AtomStereopermutator::numAssignments,
    "Returns the number of possible assignments"
  );

  atomStereopermutator.def(
    "num_stereopermutations",
    &AtomStereopermutator::numStereopermutations,
    "Returns the number of stereopermutations"
  );

  atomStereopermutator.def(pybind11::self == pybind11::self);
  atomStereopermutator.def(pybind11::self != pybind11::self);
  atomStereopermutator.def(pybind11::self < pybind11::self);

  atomStereopermutator.def(
    "__repr__",
    [](const AtomStereopermutator& permutator) {
      return permutator.info();
    }
  );
}
