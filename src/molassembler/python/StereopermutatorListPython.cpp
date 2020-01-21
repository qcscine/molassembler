/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "OptionalPython.h"

#include "molassembler/StereopermutatorList.h"

void init_stereopermutator_list(pybind11::module& m) {
  using namespace Scine::molassembler;

  pybind11::class_<StereopermutatorList> stereopermutatorList(
    m,
    "StereopermutatorList",
    R"delim(
      Manages all stereopermutators that are part of a molecule.

      >>> # A sample molecule with one stereogenic atom and bond stereopermutator each
      >>> mol = io.experimental.from_smiles("CC=C[C@](F)(Cl)[H]")
      >>> permutators = mol.stereopermutators
      >>> is_stereogenic = lambda p: p.num_assignments > 1
      >>> atom_permutators = permutators.atom_stereopermutators()
      >>> bond_permutators = permutators.bond_stereopermutators()
      >>> stereogenic_atom_permutators = [p for p in atom_permutators if is_stereogenic(p)]
      >>> stereogenic_bond_permutators = [p for p in bond_permutators if is_stereogenic(p)]
      >>> # Atom stereopermutators are instantiated on every non-terminal atom
      >>> permutators.A() > len(stereogenic_atom_permutators)
      True
      >>> len(stereogenic_atom_permutators) # But only one of them is stereogenic here
      1
      >>> # Bond stereopermutators are instantiated only where they are stereogenic
      >>> # or conserve information on relative spatial orientation
      >>> permutators.B() > len(stereogenic_bond_permutators)
      False
      >>> len(stereogenic_bond_permutators)
      1
      >>> permutators.has_unassigned_permutators() # The stereo of the double bond is unspecified
      True
      >>> permutators.has_zero_assignment_permutators()
      False
      >>> double_bond_index = BondIndex(1, 2)
      >>> assert mol.graph.bond_type(double_bond_index) == BondType.Double
      >>> # When looking up a stereopermutator, remember that you can get None back
      >>> bond_stereopermutator = permutators.option(double_bond_index)
      >>> bond_stereopermutator is None
      False
      >>> is_stereogenic(bond_stereopermutator)
      True
    )delim"
  );

  stereopermutatorList.def(
    "empty",
    &StereopermutatorList::empty,
    R"delim(
      Whether the list is empty or not. Remember that since atom
      stereopermutators are instantiated on all non-terminal atoms, this is
      rare.

      >>> h2 = Molecule() # Default constructor makes the diatomic hydrogen molecule
      >>> h2.stereopermutators.empty() # Both atoms are terminal here, so no permutators
      True
    )delim"
  );

  stereopermutatorList.def(
    "has_zero_assignment_permutators",
    &StereopermutatorList::hasZeroAssignmentStereopermutators,
    "Returns whether the list contains any stereopermutators that have no "
    "possible stereopermutations"
  );

  stereopermutatorList.def(
    "has_unassigned_permutators",
    &StereopermutatorList::hasUnassignedStereopermutators,
    "Returns whether the list contains any stereopermutators that are "
    "unassigned"
  );

  stereopermutatorList.def(
    "option",
    pybind11::overload_cast<AtomIndex>(
      &StereopermutatorList::option,
      pybind11::const_
    ),
    pybind11::arg("atom"),
    R"delim(
      Fetches a read-only option to an
      :class:`AtomStereopermutator`, if present on this atom index
    )delim"
  );

  stereopermutatorList.def(
    "option",
    pybind11::overload_cast<const BondIndex&>(
      &StereopermutatorList::option,
      pybind11::const_
    ),
    pybind11::arg("bond_index"),
    R"delim(
      Fetches a read-only option to a
      :class:`BondStereopermutator`, if present on this atom index
    )delim"
  );

  stereopermutatorList.def(
    "A",
    &StereopermutatorList::A,
    "Returns the number of :class:`AtomStereopermutator` instances stored"
  );

  stereopermutatorList.def(
    "B",
    &StereopermutatorList::B,
    "Returns the number of :class:`BondStereopermutator` instances stored"
  );

  stereopermutatorList.def(
    "atom_stereopermutators",
    [](const StereopermutatorList& list) {
      auto range = list.atomStereopermutators();
      return pybind11::make_iterator(range.begin(), range.end());
    },
    "Returns a range of all :class:`AtomStereopermutator`"
  );

  stereopermutatorList.def(
    "bond_stereopermutators",
    [](const StereopermutatorList& list) {
      auto range = list.bondStereopermutators();
      return pybind11::make_iterator(range.begin(), range.end());
    },
    "Returns a range of all :class:`BondStereopermutator`"
  );

  stereopermutatorList.def(
    "__repr__",
    [](const StereopermutatorList& list) -> std::string {
      std::string repr = "Stereopermutator list ";
      if(list.empty()) {
        repr += "(empty)";
        return repr;
      }

      const unsigned A = list.A();
      const unsigned B = list.B();

      repr += "with ";
      if(A > 0) {
        repr += std::to_string(A);
        repr += " atom ";

        if(B > 0) {
          // If B > 0, then A must be at least 2, so "and" is always justified
          repr += "and " + std::to_string(B) + " bond ";
        }

        repr += "stereopermutator";
        if(A + B > 1) {
          repr += "s";
        }
      }

      return repr;
    }
  );
}
