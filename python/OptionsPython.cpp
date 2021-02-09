/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "Molassembler/Options.h"

void init_options(pybind11::module& m) {
  using namespace Scine::Molassembler;

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
    "DoNotPreserve",
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
  pybind11::class_<Options::Thermalization> thermalization(options, "Thermalization");
  thermalization.def_readwrite_static(
    "pyramidal_inversion",
    &Options::Thermalization::pyramidalInversion,
    "When set, pyramidal nitrogen atoms not part of a small cycle invert quickly"
  );
  thermalization.def_readwrite_static(
    "berry_pseudorotation",
    &Options::Thermalization::berryPseudorotation,
    R"delim(
      If set and there are no linked substituents in a trigonal bipyramid
      shape, stereopermutations are thermalized.
    )delim"
  );
  thermalization.def_readwrite_static(
    "bartell_mechanism",
    &Options::Thermalization::bartellMechanism,
    R"delim(
      If set and there are no linked substituents in a pentagonal bipyramid
      shape, stereopermutations are thermalized.
    )delim"
  );
  thermalization.def(
    "enable",
    &Options::Thermalization::enable,
    "Sets a high temperature approximation where all thermalizations are enabled"
  );
  thermalization.def(
    "disable",
    &Options::Thermalization::disable,
    "Sets a low temperature approximation where all thermalizations are disabled"
  );

  options.def_readwrite_static(
    "chiral_state_preservation",
    &Options::chiralStatePreservation,
    "Global chiral state preservation setting of the library. Defaults to effortless and unique"
  );

  /* Access to the PRNG instance */
  m.def("randomness_engine", &randomnessEngine);
}
