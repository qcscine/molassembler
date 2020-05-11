/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

#include "Molassembler/Cycles.h"

void init_cycles(pybind11::module& m) {
  using namespace Scine::Molassembler;
  pybind11::class_<Cycles> cycles(
    m,
    "Cycles",
    R"delim(
      Information about molecular graph cycles.

      >>> # Simple molecule for which relevant cycles and cycle families are the same
      >>> spiro = io.experimental.from_smiles("C12(CCC1)CCC2")
      >>> cycles = spiro.graph.cycles
      >>> cycles.num_cycle_families()
      2
      >>> cycles.num_cycle_families(0) # The spiroatom belongs to both families
      2
      >>> cycles.num_cycle_families(1) # Other cycle atoms only belong to one
      1
      >>> cycles.num_cycle_families() == cycles.num_relevant_cycles()
      True
      >>> cycles.num_cycle_families(1) == cycles.num_relevant_cycles(1)
      True
    )delim"
  );

  cycles.def(
    "num_cycle_families",
    pybind11::overload_cast<>(
      &Cycles::numCycleFamilies,
      pybind11::const_
    ),
    "Returns the number of cycle families present in the graph"
  );
  cycles.def(
    "num_cycle_families",
    pybind11::overload_cast<AtomIndex>(
      &Cycles::numCycleFamilies,
      pybind11::const_
    ),
    pybind11::arg("constituting_index"),
    "Returns the number of cycle families an atom index belongs to"
  );
  cycles.def(
    "num_relevant_cycles",
    pybind11::overload_cast<>(
      &Cycles::numRelevantCycles,
      pybind11::const_
    ),
    "Returns the number of relevant cycles present in the graph"
  );
  cycles.def(
    "num_relevant_cycles",
    pybind11::overload_cast<AtomIndex>(
      &Cycles::numRelevantCycles,
      pybind11::const_
    ),
    pybind11::arg("constituting_index"),
    "Returns the number of relevant cycles an atom index belongs to"
  );

  cycles.def(
    "__iter__",
    [](const Cycles& cyc) {return pybind11::make_iterator(cyc.begin(), cyc.end());},
    "Iterate through all relevant cycles."
  );

  cycles.def(
    "__len__",
    pybind11::overload_cast<>(
      &Cycles::numRelevantCycles,
      pybind11::const_
    ),
    "Returns the number of relevant cycles in the graph"
  );
}
