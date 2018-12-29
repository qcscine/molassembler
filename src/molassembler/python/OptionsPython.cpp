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
  ).value("Low", TemperatureRegime::Low)
    .value("High", TemperatureRegime::High);

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
  ).value("None", ChiralStatePreservation::None)
    .value("EffortlessAndUnique", ChiralStatePreservation::EffortlessAndUnique)
    .value("Unique", ChiralStatePreservation::Unique)
    .value("RandomFromMultipleBest", ChiralStatePreservation::RandomFromMultipleBest);

  /* Tau criterion */
  pybind11::enum_<TauCriterion>(
    m,
    "TauCriterion",
    "Specifies use of the tau criterion in differentiating between symmetries "
    "of sizes four and five."
  ).value("Enable", TauCriterion::Enable)
    .value("Disable", TauCriterion::Disable);

  pybind11::class_<Options> options(m, "Options", "Contains global library settings");
  options.def_readwrite_static("temperature_regime", &Options::temperatureRegime);
  options.def_readwrite_static("chiral_state_preservation", &Options::chiralStatePreservation);
  options.def_readwrite_static("tau_criterion", &Options::tauCriterion);

  /* Access to the PRNG instance */
  m.def("randomness_engine", &randomnessEngine);
}
