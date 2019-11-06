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
    R"delim(
      This class handles the permutation of ranked ligands around a central
      atom. It models its haptic ligands' binding sites and bridges in
      multidentate ligands in order to decide which stereopermutations are
      feasible. A stereopermutation may be infeasible, i.e. not realizable in
      three-dimensional space, if either haptic ligands would intersect due to
      too close ligand angles for their spatial heft, or if a multidentate
      ligand's bridge length between binding sites were too short to match the
      angle. The list of stereopermutations reduced by infeasible
      stereopermutations is then re-indexed and those indices into the list are
      called assignments.

      A Stereopermutator can be unassigned, i.e. the distinct stereopermutation
      that the substituents are can be indeterminate. If you choose to generate
      conformers for a molecule that includes unassigned stereopermutators,
      every conformer will choose an assignment from the pool of feasible
      assignments randomly, but consistent with relative statistical occurrence
      weights.
    )delim"
  );

  atomStereopermutator.def(
    "angle",
    &AtomStereopermutator::angle,
    pybind11::arg("site_index_i"),
    pybind11::arg("site_index_j"),
    "Fetches the angle between substituent site indices in radians"
  );

  atomStereopermutator.def(
    "assigned",
    &AtomStereopermutator::assigned,
    "Returns the assignment integer if assigned, ``None`` otherwise."
  );

  atomStereopermutator.def(
    "central_index",
    &AtomStereopermutator::centralIndex,
    "Returns the central atom this permutator is placed on"
  );

  atomStereopermutator.def(
    "index_of_permutation",
    &AtomStereopermutator::indexOfPermutation,
    R"delim(
      Returns the index of permutation if assigned, otherwise returns ``None``.
    )delim"
  );

  atomStereopermutator.def_property_readonly(
    "ranking",
    &AtomStereopermutator::getRanking,
    R"delim(
      Get the underlying ranking state of substituents

      :rtype: :class:`RankingInformation`
    )delim"
  );

  atomStereopermutator.def_property_readonly(
    "shape",
    &AtomStereopermutator::getShape,
    R"delim(
      Returns the underlying shape

      :rtype: :class:`shapes.Shape`
    )delim"
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
