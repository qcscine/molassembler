/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "molassembler/Cycles.h"

void init_cycles(pybind11::module& m) {
  using namespace Scine::molassembler;
  pybind11::class_<Cycles> cycles(m, "Cycles", "Information about molecular graph cycles.");

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
