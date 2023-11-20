/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/Descriptors.h"

void init_descriptors(pybind11::module& m) {
  using namespace Scine::Molassembler;

  m.def(
    "ranking_equivalent_groups",
    &rankingEquivalentGroups,
    pybind11::arg("molecule"),
    R"delim(
      Determines ranking equivalent groups of atoms

      Uses atom-stereopermutator ranking information to determine which parts of
      molecules are completely ranking-equivalent. Note that ranking symmetry at
      bonds is not found by this method, yielding more ranking equivalent atoms
      than there truly are in many cases, such as in ethane.

      :returns: A list of length equal to the number of atoms in the molecule
        containing group indices. Atoms with the same group index are ranking
        equivalent.

      >>> cyclopentane = io.experimental.from_smiles("C1CCCC1")
      >>> groups = ranking_equivalent_groups(cyclopentane)
      >>> groups
      [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
      >>> len(groups) == cyclopentane.graph.N
      True
    )delim"
  );

  m.def(
    "ranking_distinct_atoms",
    &rankingDistinctAtoms,
    pybind11::arg("molecule"),
    R"delim(
      Determines the set of atoms not ranked equivalently at a permutator

      Uses atom-stereopermutator ranking information to determine which parts of
      molecules are completely ranking-equivalent. Note that ranking symmetry at
      bonds is not found by this method, yielding more ranking equivalent atoms
      than there truly are in many cases, such as in ethane.

      :returns: An unordered list of atoms

      >>> propane = io.experimental.from_smiles("CCC")
      >>> propane.graph.N
      11
      >>> distinct = sorted(ranking_distinct_atoms(propane))
      >>> distinct
      [0, 1, 3, 6]
    )delim"
  );
}
