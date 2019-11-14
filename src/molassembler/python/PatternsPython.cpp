/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"

#include "molassembler/Patterns.h"
#include "molassembler/Molecule.h"

void init_patterns_submodule(pybind11::module& m) {
  using namespace Scine;
  using namespace molassembler;

  auto patternsSubmodule = m.def_submodule("patterns");
  patternsSubmodule.doc() = R"(Patterns submodule)";

  patternsSubmodule.def(
    "alkane",
    &patterns::alkane,
    pybind11::arg("N"),
    R"delim(
      Generate an alkane of length N

      >>> alkane_atom_count = lambda n: 4 if n == 1 else 3 * n + 2
      >>> methane = alkane(1)
      >>> methane.graph.N == alkane_atom_count(1)
      True
      >>> alkane(4).graph.N == alkane_atom_count(4)
      True
    )delim"
  );
}
