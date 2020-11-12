/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/Descriptors.h"

void init_descriptors(pybind11::module& m) {
  using namespace Scine::Molassembler;

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
    )delim"
  );
}
