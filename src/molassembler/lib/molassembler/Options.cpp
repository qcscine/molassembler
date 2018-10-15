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

bool disregardStereopermutator(
  const AtomStereopermutator& stereopermutator,
  const Delib::ElementType centralType,
  const Cycles& cycleData,
  const TemperatureRegime temperatureRegimeSetting
) {
  if(
    temperatureRegimeSetting == TemperatureRegime::High
    && stereopermutator.getSymmetry() == Symmetry::Name::CutTetrahedral
    && centralType == Delib::ElementType::N
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

void pickyFit(
  AtomStereopermutator& stereopermutator,
  const OuterGraph& graph,
  const AngstromWrapper& angstromWrapper,
  const Symmetry::Name expectedSymmetry
) {
  /* Seesaw, trigonal pyramidal and tetrahedral are surprisingly close in terms
   * of angles, and sometimes just slightly distorted tetrahedral centers can
   * be recognized as seesaws, even though it makes absolutely zero sense. So
   * in case the atom is a carbon, the expected geometry is tetrahedral and it
   * has four adjacencies, just exclude Seesaw and trigonal pyramidal from the
   * list of symmetries being fitted against.
   *
   * Calling
   * determineLocalGeometry is somewhat overkill here, but possibly more
   * future-proof.
   */
  if(
    graph.elementType(stereopermutator.centralIndex()) == Delib::ElementType::C
    && expectedSymmetry == Symmetry::Name::Tetrahedral
  ) {
    stereopermutator.fit(
      graph,
      angstromWrapper,
      {Symmetry::Name::Seesaw, Symmetry::Name::TrigonalPyramidal}
    );
  } else {
    stereopermutator.fit(graph, angstromWrapper);
  }
}


} // namespace molassembler
