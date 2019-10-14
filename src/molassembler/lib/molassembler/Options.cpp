/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Options.h"

#include "boost/range/iterator_range_core.hpp"
#include "chemical_symmetries/Symmetries.h"

#include "molassembler/AtomStereopermutator.h"
#include "molassembler/Cycles.h"
#include "molassembler/OuterGraph.h"

namespace Scine {

namespace molassembler {

random::Engine& randomnessEngine() {
  // Pursuant to Construct-on-first-use idiom
  static random::Engine engine;
  return engine;
}

TemperatureRegime Options::temperatureRegime = TemperatureRegime::High;
ChiralStatePreservation Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;
TauCriterion Options::tauCriterion = TauCriterion::Enable;
ShapeTransition Options::shapeTransition = ShapeTransition::MaximizeChiralStatePreservation;

bool disregardStereopermutator(
  const AtomStereopermutator& stereopermutator,
  const Scine::Utils::ElementType centralType,
  const Cycles& cycleData,
  const TemperatureRegime temperatureRegimeSetting
) {
  if(
    temperatureRegimeSetting == TemperatureRegime::High
    && stereopermutator.getShape() == Symmetry::Shape::ApicalTrigonalPyramid
    && centralType == Scine::Utils::ElementType::N
  ) {
    // Figure out if the nitrogen is in a cycle of size 4 or smaller
    for(
      const auto cycleEdges :
      boost::make_iterator_range(
        cycleData.containing(
          stereopermutator.centralIndex()
        )
      )
    ) {
      if(cycleEdges.size() <= 4) {
        return false;
      }
    }

    return true;
  }

  return false;
}

} // namespace molassembler

} // namespace Scine
