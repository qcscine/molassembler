/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "molassembler/Options.h"

void init_options(pybind11::module& m) {
  using namespace Scine::molassembler;

  /* Temperature regime */
  pybind11::enum_<TemperatureRegime>(
    m,
    "TemperatureRegime",
    "Temperature we model molecules at. Low means that no pyramidal inversion "
    "occurs, so that all nitrogen atoms can be stereocenters. If High is set, "
    "then only nitrogen geometries in particularly strained cycles (3, 4) can "
    "be stereocenters"
  ).value("Low", TemperatureRegime::Low, "No pyramidal inversion, all nitrogen atoms can be stereopermutators")
    .value("High", TemperatureRegime::High, "Only nitrogen atoms in particularly strained cycles (3, 4) can be stereopermutators");

  /* Chiral state preservation */
  pybind11::enum_<ChiralStatePreservation>(
    m,
    "ChiralStatePreservation",
    "Specifies how chiral state is to be propagated on graph modifications.  If "
    "None is set, no chiral state is preserved. If EffortlessAndUnique is set, "
    "only unambiguous zero-effort mappings are used to propagate chiral state. "
    "This enables e.g. the propagation of ligand loss in octahedral to square "
    "pyramidal and back.  If Unique is set, chiral state is propagated if the "
    "best mapping is unique, i.e. there are no other mappings with the same "
    "quality measures. Enables e.g. the propagation of seesaw to square planar, "
    "but not back. Under RandomFromMultipleBest, random mappings are chosen "
    "from the set of best mappings, permitting chiral state propagation in all "
    "cases. "
  ).value("None", ChiralStatePreservation::None, "Don't try to preserve chiral state")
    .value(
      "EffortlessAndUnique",
      ChiralStatePreservation::EffortlessAndUnique,
      "Use only completely unambiguous zero-effort mappings. Those are the green"
      "edges in the graphs below. Note that the ligand gain situation from square"
      "planar to square pyramidal is not unique, and therefore not shown as"
      "green."
    ).value(
      "Unique",
      ChiralStatePreservation::Unique,
      "Propagates if the best shape mapping is unique, i.e. there are no other"
      "mappings with the same quality measures. This enables all green and black"
      "edges."
    ).value(
      "RandomFromMultipleBest",
      ChiralStatePreservation::RandomFromMultipleBest,
      "Chooses randomly from the set of best mappings, permitting chiral state"
      "propagation in all cases. So propagating chiral state from square planar"
      "to square pyramidal is now possible -- there are two ways of placing the"
      "new apical ligand -- but you only get one of them."
    );

  pybind11::class_<Options> options(m, "Options", "Contains global library settings");
  options.def_readwrite_static(
    "temperature_regime",
    &Options::temperatureRegime,
    "Global temperature regime setting of the library. Defaults to high temperature approximation"
  );
  options.def_readwrite_static(
    "chiral_state_preservation",
    &Options::chiralStatePreservation,
    "Global chiral state preservation setting of the library. Defaults to effortless and unique"
  );

  /* Access to the PRNG instance */
  m.def("randomness_engine", &randomnessEngine);
}
