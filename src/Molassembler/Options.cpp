/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Options.h"

#include "Molassembler/Shapes/Data.h"

namespace Scine {
namespace Molassembler {

Random::Engine& randomnessEngine() {
  // Pursuant to Construct-on-first-use idiom
  static Random::Engine engine;
  return engine;
}

bool Options::Thermalization::pyramidalInversion = true;
bool Options::Thermalization::berryPseudorotation = true;
bool Options::Thermalization::bartellMechanism = true;

ChiralStatePreservation Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;
ShapeTransition Options::shapeTransition = ShapeTransition::MaximizeChiralStatePreservation;

} // namespace Molassembler
} // namespace Scine
