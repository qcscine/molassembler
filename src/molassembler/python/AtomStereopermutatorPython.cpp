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

      Stereopermutator instances themselves are nonmodifiable. To change
      them, you have to make changes at the molecule level.

      >>> import molassembler as masm
      >>> methane = masm.patterns.alkane(1)
      >>> methane_central_stereopermutator = methane.stereopermutators.option(0)
      >>> methane_central_stereopermutator is not None
      True
      >>> methane_central_stereopermutator.shape == masm.shapes.Shape.Tetrahedron
      True
      >>> methane_central_stereopermutator
      A on 0 (tetrahedron, AAAA): 0/1
    )delim"
  );

  atomStereopermutator.def(
    "angle",
    &AtomStereopermutator::angle,
    pybind11::arg("site_index_i"),
    pybind11::arg("site_index_j"),
    R"delim(
      Fetches the angle between substituent site indices in radians

      >>> import math
      >>> tetrahedron_angle = 2 * math.atan(math.sqrt(2))
      >>> import molassembler as masm
      >>> methane = masm.patterns.alkane(1)
      >>> a = methane.stereopermutators.option(0)
      >>> math.isclose(a.angle(0, 1), tetrahedron_angle)
      True
    )delim"
  );

  atomStereopermutator.def_property_readonly(
    "assigned",
    &AtomStereopermutator::assigned,
    "The assignment integer if assigned, ``None`` otherwise."
  );

  atomStereopermutator.def_property_readonly(
    "central_index",
    &AtomStereopermutator::centralIndex,
    "The central atom this permutator is placed on"
  );

  atomStereopermutator.def_property_readonly(
    "index_of_permutation",
    &AtomStereopermutator::indexOfPermutation,
    R"delim(
      The index of permutation if assigned, otherwise ``None``.
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

  atomStereopermutator.def_property_readonly(
    "num_assignments",
    &AtomStereopermutator::numAssignments,
    "The number of feasible assignments"
  );

  atomStereopermutator.def_property_readonly(
    "num_stereopermutations",
    &AtomStereopermutator::numStereopermutations,
    "The number of stereopermutations"
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
