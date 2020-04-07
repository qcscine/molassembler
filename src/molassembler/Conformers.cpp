/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Temple/Functional.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

namespace Scine {

namespace molassembler {

std::vector<
  outcome::result<Utils::PositionCollection>
> generateRandomEnsemble(
  const Molecule& molecule,
  const unsigned numStructures,
  const distance_geometry::Configuration& configuration
) {
  auto result = distance_geometry::run(molecule, numStructures, configuration, boost::none);

  /* Convert the AngstromPositionss into PositionCollections */
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

std::vector<
  outcome::result<Utils::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  const unsigned numStructures,
  const unsigned seed,
  const distance_geometry::Configuration& configuration
) {
  auto result = distance_geometry::run(molecule, numStructures, configuration, seed);

  /* Convert the AngstromPositionss into PositionCollections */
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

outcome::result<Utils::PositionCollection> generateRandomConformation(
  const Molecule& molecule,
  const distance_geometry::Configuration& configuration
) {
  auto result = distance_geometry::run(molecule, 1, configuration, boost::none);

  assert(result.size() == 1);
  auto& wrapperResult = result.front();

  if(wrapperResult) {
    return wrapperResult.value().getBohr();
  }

  return wrapperResult.as_failure();
}

outcome::result<Utils::PositionCollection> generateConformation(
  const Molecule& molecule,
  const unsigned seed,
  const distance_geometry::Configuration& configuration
) {
  auto result = distance_geometry::run(molecule, 1, configuration, seed);

  assert(result.size() == 1);
  auto& wrapperResult = result.front();

  if(wrapperResult) {
    return wrapperResult.value().getBohr();
  }

  return wrapperResult.as_failure();
}

} // namespace molassembler

} // namespace Scine
