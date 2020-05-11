/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Options.h"

#include "molassembler/Shapes/Data.h"

namespace Scine {
namespace Molassembler {

Random::Engine& randomnessEngine() {
  // Pursuant to Construct-on-first-use idiom
  static Random::Engine engine;
  return engine;
}

TemperatureRegime Options::temperatureRegime = TemperatureRegime::High;
ChiralStatePreservation Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;
ShapeTransition Options::shapeTransition = ShapeTransition::MaximizeChiralStatePreservation;

} // namespace Molassembler
} // namespace Scine
