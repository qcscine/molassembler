/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Options.h"

#include "boost/range/iterator_range_core.hpp"
#include "shapes/Data.h"

namespace Scine {

namespace molassembler {

random::Engine& randomnessEngine() {
  // Pursuant to Construct-on-first-use idiom
  static random::Engine engine;
  return engine;
}

TemperatureRegime Options::temperatureRegime = TemperatureRegime::High;
ChiralStatePreservation Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;
ShapeTransition Options::shapeTransition = ShapeTransition::MaximizeChiralStatePreservation;

} // namespace molassembler

} // namespace Scine
