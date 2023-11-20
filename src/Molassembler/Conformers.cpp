/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Temple/Functional.h"
#include "Molassembler/DistanceGeometry/ConformerGeneration.h"

namespace Scine {
namespace Molassembler {

std::vector<
  Result<Utils::PositionCollection>
> generateRandomEnsemble(
  const Molecule& molecule,
  const unsigned numStructures,
  const DistanceGeometry::Configuration& configuration
) {
  auto result = DistanceGeometry::run(molecule, numStructures, configuration, boost::none);

  /* Convert the AngstromPositionss into PositionCollections */
  std::vector<
    Result<Utils::PositionCollection>
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

std::vector<
  Result<Utils::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  const unsigned numStructures,
  const unsigned seed,
  const DistanceGeometry::Configuration& configuration
) {
  auto result = DistanceGeometry::run(molecule, numStructures, configuration, seed);

  /* Convert the AngstromPositionss into PositionCollections */
  std::vector<
    Result<Utils::PositionCollection>
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

Result<Utils::PositionCollection> generateRandomConformation(
  const Molecule& molecule,
  const DistanceGeometry::Configuration& configuration
) {
  auto result = DistanceGeometry::run(molecule, 1, configuration, boost::none);

  assert(result.size() == 1);
  auto& wrapperResult = result.front();

  if(wrapperResult) {
    return wrapperResult.value().getBohr();
  }

  return wrapperResult.as_failure();
}

Result<Utils::PositionCollection> generateConformation(
  const Molecule& molecule,
  const unsigned seed,
  const DistanceGeometry::Configuration& configuration
) {
  auto result = DistanceGeometry::run(molecule, 1, configuration, seed);

  assert(result.size() == 1);
  auto& wrapperResult = result.front();

  if(wrapperResult) {
    return wrapperResult.value().getBohr();
  }

  return wrapperResult.as_failure();
}

} // namespace Molassembler
} // namespace Scine
