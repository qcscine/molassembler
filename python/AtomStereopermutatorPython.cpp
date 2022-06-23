/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "Molassembler/Temple/Functional.h"
#include "TypeCasters.h"
#include "pybind11/operators.h"

#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/Graph.h"

void init_atom_stereopermutator(pybind11::module& m) {
  using namespace Scine::Molassembler;
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

      >>> # Enantiomeric pair of asymmetric tetrahedral carbon atoms
      >>> import scine_utilities as utils
      >>> asym_carbon = io.experimental.from_smiles("N[C@](Br)(O)F")
      >>> carbon_index = asym_carbon.graph.atoms_of_element(utils.ElementType.C)[0]
      >>> carbon_stereopermutator = asym_carbon.stereopermutators.option(carbon_index)
      >>> assert carbon_stereopermutator is not None
      >>> carbon_stereopermutator.shape == shapes.Shape.Tetrahedron
      True
      >>> carbon_stereopermutator.assigned is not None
      True
      >>> enantiomer = io.experimental.from_smiles("N[C@@](Br)(O)F")
      >>> assert enantiomer.graph.element_type(carbon_index) == utils.ElementType.C
      >>> enantiomer_stereopermutator = enantiomer.stereopermutators.option(carbon_index)
      >>> enantiomer_stereopermutator.assigned is not None
      True
      >>> carbon_stereopermutator.assigned != enantiomer_stereopermutator.assigned
      True
    )delim"
  );

  atomStereopermutator.def(
    pybind11::init<const Graph&, Shapes::Shape, AtomIndex, RankingInformation>(),
    pybind11::arg("graph"),
    pybind11::arg("shape"),
    pybind11::arg("placement"),
    pybind11::arg("ranking"),
    "Instantiate an atom stereopermutator at a particular position in a graph"
  );

  atomStereopermutator.def(
    "angle",
    &AtomStereopermutator::angle,
    pybind11::arg("site_index_i"),
    pybind11::arg("site_index_j"),
    R"delim(
      Fetches the angle between substituent site indices in radians

      >>> # The tetrahedron angle
      >>> import math
      >>> tetrahedron_angle = 2 * math.atan(math.sqrt(2))
      >>> methane = io.experimental.from_smiles("C")
      >>> a = methane.stereopermutators.option(0)
      >>> math.isclose(a.angle(0, 1), tetrahedron_angle)
      True
    )delim"
  );

  atomStereopermutator.def_property_readonly(
    "assigned",
    &AtomStereopermutator::assigned,
    R"delim(
      The assignment integer if assigned, ``None`` otherwise.

      >>> # A stereo-unspecified tetrahedral asymmetric carbon atom
      >>> import scine_utilities as utils
      >>> asymmetric_carbon = io.experimental.from_smiles("NC(Br)(O)F")
      >>> carbon_index = asymmetric_carbon.graph.atoms_of_element(utils.ElementType.C)[0]
      >>> stereopermutator = asymmetric_carbon.stereopermutators.option(carbon_index)
      >>> assert stereopermutator is not None
      >>> stereopermutator.num_assignments == 2
      True
      >>> stereopermutator.assigned is None
      True
    )delim"
  );

  atomStereopermutator.def_property_readonly(
    "placement",
    &AtomStereopermutator::placement,
    "The central atom this permutator is placed on"
  );

  atomStereopermutator.def_property_readonly(
    "index_of_permutation",
    &AtomStereopermutator::indexOfPermutation,
    R"delim(
      The index of permutation if assigned, otherwise ``None``. Indices
      of permutation are the abstract index of permutation within the set of
      permutations that do not consider feasibility. This is not necessarily
      equal to the assignment index.

      >>> # The shipscrew lambda and delta isomers where trans-ligation is impossible
      >>> shipscrew_smiles = "[Fe@OH1+3]123(OC(=O)C(=O)O1)(OC(=O)C(=O)O2)OC(=O)C(=O)O3"
      >>> shipscrew = io.experimental.from_smiles(shipscrew_smiles)
      >>> permutator = shipscrew.stereopermutators.option(0)
      >>> assert permutator is not None
      >>> permutator.num_stereopermutations # Number of abstract permutations
      4
      >>> permutator.num_assignments # Number of spatially feasible permutations
      2
      >>> permutator.index_of_permutation
      2
      >>> permutator.assigned
      1
      >>> shipscrew.assign_stereopermutator(0, None) # Dis-assign the stereopermutator
      >>> permutator = shipscrew.stereopermutators.option(0)
      >>> assert permutator is not None
      >>> permutator.index_of_permutation is None
      True
      >>> permutator.assigned is None
      True
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
    "thermalized",
    pybind11::overload_cast<AtomIndex, Shapes::Shape, const RankingInformation&, const Graph&>(&AtomStereopermutator::thermalized),
    "Whether the stereopermutations are thermalized"
  );

  atomStereopermutator.def_property_readonly(
    "num_assignments",
    &AtomStereopermutator::numAssignments,
    "The number of feasible assignments. See index_of_permutation."
  );

  atomStereopermutator.def_property_readonly(
    "num_stereopermutations",
    &AtomStereopermutator::numStereopermutations,
    "The number of stereopermutations. See index_of_permutation."
  );

  atomStereopermutator.def(
    "site_groups",
    &AtomStereopermutator::siteGroups,
    R"delim(
      Site indices grouped by whether their shape vertices are interconvertible
      by rotation. Can be used to distinguish e.g. apical / equatorial ligands
      in a square pyramid shape.

      :raises: RuntimeError if the stereopermutator is not assigned
    )delim"
  );

  atomStereopermutator.def_property_readonly(
    "vertex_map",
    [](const AtomStereopermutator& perm) -> std::vector<Shapes::Vertex> {
      return Temple::map(perm.getShapePositionMap(), Temple::Functor::second);
    }
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
