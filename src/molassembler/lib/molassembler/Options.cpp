#include "Options.h"

#include "CNStereocenter.h"
#include "Cycles.h"
#include "GraphHelpers.h"

namespace molassembler {

TemperatureRegime Options::temperatureRegime = TemperatureRegime::High;
ChiralStatePreservation Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;

bool disregardStereocenter(
  const Stereocenters::CNStereocenter& stereocenter,
  const Delib::ElementType centralType,
  const Cycles& cycleData,
  const TemperatureRegime temperatureRegimeSetting
) {
  if(temperatureRegimeSetting == TemperatureRegime::High) {
    AtomIndexType centralIndex = stereocenter.involvedAtoms().front();

    if(
      stereocenter.getSymmetry() == Symmetry::Name::CutTetrahedral
      && centralType == Delib::ElementType::N
    ) {
      // Figure out if the nitrogen is in a cycle of size 4 or smaller
      for(
        const auto cyclePtr :
        cycleData.iterate(Cycles::predicates::ContainsIndex {centralIndex})
      ) {
        if(Cycles::size(cyclePtr) <= 4) {
          return false;
        }
      }

      return true;
    }

    return false;
  }

  return false;
}

void pickyFit(
  Stereocenters::CNStereocenter& stereocenter,
  const GraphType& graph,
  const AngstromWrapper& angstromWrapper,
  const Symmetry::Name expectedSymmetry
) {
  const AtomIndexType centralIndex = stereocenter.involvedAtoms().front();

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
    graph::elementType(centralIndex, graph) == Delib::ElementType::C
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
