// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Options.h"

#include "boost/range/iterator_range_core.hpp"
#include "chemical_symmetries/Symmetries.h"

#include "molassembler/AtomStereopermutator.h"
#include "molassembler/Cycles.h"
#include "molassembler/OuterGraph.h"

namespace molassembler {

random::Engine& randomnessEngine() {
  // Pursuant to Construct-on-first-use idiom
  static random::Engine engine;
  return engine;
}

TemperatureRegime Options::temperatureRegime = TemperatureRegime::High;
ChiralStatePreservation Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;
TauCriterion Options::tauCriterion = TauCriterion::Enable;

bool disregardStereopermutator(
  const AtomStereopermutator& stereopermutator,
  const Scine::Utils::ElementType centralType,
  const Cycles& cycleData,
  const TemperatureRegime temperatureRegimeSetting
) {
  if(
    temperatureRegimeSetting == TemperatureRegime::High
    && stereopermutator.getSymmetry() == Symmetry::Name::CutTetrahedral
    && centralType == Scine::Utils::ElementType::N
  ) {
    // Figure out if the nitrogen is in a cycle of size 4 or smaller
    for(
      const auto cyclePtr :
      boost::make_iterator_range(
        cycleData.iteratorPair(
          Cycles::predicates::ContainsIndex {stereopermutator.centralIndex()}
        )
      )
    ) {
      if(Cycles::size(cyclePtr) <= 4) {
        return false;
      }
    }

    return true;
  }

  return false;
}

} // namespace molassembler
