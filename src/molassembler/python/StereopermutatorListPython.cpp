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
    "Manages all stereopermutators that are part of a molecule"
  );

  stereopermutatorList.def(
    "empty",
    &StereopermutatorList::empty,
    "Whether the list is empty or not"
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
    "Returns the number of :class:`AtomStereopermutator`"
  );

  stereopermutatorList.def(
    "B",
    &StereopermutatorList::B,
    "Returns the number of :class:`BondStereopermutator`"
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
}
