#include "molassembler/Options.h"

#include "chemical_symmetries/Symmetries.h"

#include "molassembler/AtomStereocenter.h"
#include "molassembler/Cycles.h"
#include "molassembler/OuterGraph.h"

namespace molassembler {

temple::Generator prng;

TemperatureRegime Options::temperatureRegime = TemperatureRegime::High;
ChiralStatePreservation Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;

bool disregardStereocenter(
  const AtomStereocenter& stereocenter,
  const Delib::ElementType centralType,
  const Cycles& cycleData,
  const TemperatureRegime temperatureRegimeSetting
) {
  if(
    temperatureRegimeSetting == TemperatureRegime::High
    && stereocenter.getSymmetry() == Symmetry::Name::CutTetrahedral
    && centralType == Delib::ElementType::N
  ) {
    // Figure out if the nitrogen is in a cycle of size 4 or smaller
    for(
      const auto cyclePtr :
      cycleData.iterate(Cycles::predicates::ContainsIndex {stereocenter.centralIndex()})
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
  AtomStereocenter& stereocenter,
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
    graph.elementType(stereocenter.centralIndex()) == Delib::ElementType::C
    && expectedSymmetry == Symmetry::Name::Tetrahedral
  ) {
    stereocenter.fit(
      graph,
      angstromWrapper,
      {Symmetry::Name::Seesaw, Symmetry::Name::TrigonalPyramidal}
    );
  } else {
    stereocenter.fit(graph, angstromWrapper);
  }
}


} // namespace molassembler
