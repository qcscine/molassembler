// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Conformers.h"

#include "temple/Functional.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

namespace molassembler {

outcome::result<
  std::vector<Delib::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  const unsigned numStructures
) {
  auto result = DistanceGeometry::run(molecule, numStructures);

  if(result) {
    return temple::map(
      result.value(),
      [](AngstromWrapper wrapper) -> Delib::PositionCollection {
        return wrapper.getBohr();
      }
    );
  }

  return result.as_failure();
}

outcome::result<Delib::PositionCollection> generateConformation(const Molecule& molecule) {
  auto result = DistanceGeometry::run(molecule, 1);

  if(result) {
    auto& conformationList = result.value();
    assert(conformationList.size() == 1);
    auto& wrapper = conformationList.front();

    return wrapper.getBohr();
  }

  return result.as_failure();
}

} // namespace molassembler
