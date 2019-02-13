/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "temple/Functional.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

namespace Scine {

namespace molassembler {

std::vector<
  outcome::result<Utils::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  const unsigned numStructures,
  const DistanceGeometry::Configuration& configuration
) {
  auto result = DistanceGeometry::run(molecule, numStructures, configuration);

  /* Convert the AngstromWrappers into PositionCollections */
  std::vector<
    outcome::result<Utils::PositionCollection>
  > converted;
  converted.reserve(numStructures);

  for(auto& positionResult : result) {
    if(positionResult) {
      converted.emplace_back(
        positionResult.value().getBohr()
      );
    } else {
      converted.emplace_back(positionResult.as_failure());
    }
  }

  return converted;
}

outcome::result<Utils::PositionCollection> generateConformation(
  const Molecule& molecule,
  const DistanceGeometry::Configuration& configuration
) {
  auto result = DistanceGeometry::run(molecule, 1, configuration);

  assert(result.size() == 1);
  auto& wrapperResult = result.front();

  if(wrapperResult) {
    return wrapperResult.value().getBohr();
  }

  return wrapperResult.as_failure();
}

} // namespace molassembler

} // namespace Scine
