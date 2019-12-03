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
    R"delim(
      Temperature we model molecules at. Low means that no stereopermutations
      interconvert thermally. There is no pyramidal inversion and all nitrogen
      atoms can be stereocenters. No Berry pseudorotations or Bartell mechanisms
      occur, so trigonal bipyramid and pentagonal bipyramid centers can be
      stereocenters. If ``High`` is set, then nitrogen atoms in
      particularly strained cycles (3, 4) can be stereocenters. Berry
      pseudorotations and Bartell mechanisms thermalize all stereopermutations
      of their respective shapes if none of the substituents are linked.
    )delim"
  ).value("Low", TemperatureRegime::Low, "No stereopermutations interconvert thermally.")
    .value("High", TemperatureRegime::High, "Under specific circumstances, stereopermutations interconvert rapidly.");

  /* Chiral state preservation */
  pybind11::enum_<ChiralStatePreservation>(
    m,
    "ChiralStatePreservation",
    R"delim(
      Specifies how chiral state is to be propagated on graph modifications. If
      ``None`` is set, no chiral state is preserved. If ``EffortlessAndUnique``
      is set, only unambiguous zero-effort mappings are used to propagate
      chiral state. This enables e.g. the propagation of ligand loss in
      octahedral to square pyramidal and back.  If ``Unique`` is set, chiral
      state is propagated if the best mapping is unique, i.e. there are no
      other mappings with the same quality measures. Enables e.g. the
      propagation of seesaw to square planar, but not back. Under
      ``RandomFromMultipleBest``, random mappings are chosen from the set of
      best mappings, permitting chiral state propagation in all cases.
    )delim"
  ).value(
    "None",
    ChiralStatePreservation::None,
    R"delim(
      Don't try to preserve chiral state. Changes at stereopermutators always
      result in loss of chiral state.
    )delim"
  ).value(
    "EffortlessAndUnique",
    ChiralStatePreservation::EffortlessAndUnique,
    R"delim(
      Use only completely unambiguous zero-effort mappings. Note that for
      instance the ligand gain situation from square planar to square pyramidal
      is not unique, and therefore chiral state is not propagated there under
      this option.
    )delim"
  ).value(
    "Unique",
    ChiralStatePreservation::Unique,
    R"delim(
      Propagates if the best shape mapping is unique, i.e. there are no other
      mappings with the same quality measures.
    )delim"
  ).value(
    "RandomFromMultipleBest",
    ChiralStatePreservation::RandomFromMultipleBest,
    R"delim(
      Chooses randomly from the set of best mappings, permitting chiral state
      propagation in all cases. So propagating chiral state from square planar
      to square pyramidal is now possible -- there are two ways of placing the
      new apical ligand -- but you only get one of them.
    )delim"
  );

  pybind11::class_<Options> options(m, "Options", "Contains global library settings");
  options.def_readwrite_static(
    "temperature_regime",
    &Options::temperatureRegime,
    "Global temperature regime setting of the library. Defaults to high."
  );
  options.def_readwrite_static(
    "chiral_state_preservation",
    &Options::chiralStatePreservation,
    "Global chiral state preservation setting of the library. Defaults to effortless and unique"
  );

  /* Access to the PRNG instance */
  m.def("randomness_engine", &randomnessEngine);
}
