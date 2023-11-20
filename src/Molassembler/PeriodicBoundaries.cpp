/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/PeriodicBoundaries.h"

#include "Molassembler/Graph.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Adaptors/SequentialPairs.h"

#include <cassert>

namespace Scine {
namespace Molassembler {

SubstitutionsGenerator::SubstitutionMap
SubstitutionsGenerator::removeGhosts(
  Graph& graph,
  const PeriodicBoundaryDuplicates& periodics
) {
  SubstitutionMap replacements;
  const auto insertOrAppend = [&replacements](AtomIndex v, AtomIndex r, AtomIndex g) {
    const auto findIter = replacements.find(v);
    if(findIter == std::end(replacements)) {
      replacements.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(v),
        std::forward_as_tuple(1, std::make_pair(r, g))
      );
    } else {
      findIter->second.emplace_back(r, g);
    }
  };

  for(const auto& ghostRealPair : periodics.ghostAtomMap) {
    const AtomIndex ghost = ghostRealPair.first;
    const AtomIndex real = ghostRealPair.second;
    for(const AtomIndex ghostAdjacent : graph.adjacents(ghost)) {
      if(periodics.uninterestingAtoms.count(ghostAdjacent) == 0) {
        insertOrAppend(ghostAdjacent, real, ghost);
      }
    }
  }

  // Check preconditions prior to ghost atom removals from the graph
  auto ghosts = Temple::map(
    periodics.ghostAtomMap,
    [](const auto& ghostRealPair) -> AtomIndex {
      return ghostRealPair.first;
    }
  );

  std::sort(std::begin(ghosts), std::end(ghosts), std::greater<>());
  const bool sequentialGroup = Temple::all_of(
    Temple::Adaptors::sequentialPairs(ghosts),
    [](const AtomIndex i, const AtomIndex j) -> bool {
      return i == j + 1;
    }
  );
  const bool largestGhostIsV = (
    ghosts.empty()
    || ghosts.front() == graph.V() - 1
  );
  if(!sequentialGroup || !largestGhostIsV) {
    throw std::runtime_error("Violated pbc precondition that ghost atoms must be a contiguous block of atom indices larger than all others");
  }

  for(const AtomIndex ghost : ghosts) {
    graph.inner().clearVertex(ghost);
    graph.inner().removeVertex(ghost);
  }

  return replacements;
}

} // namespace Molassembler
} // namespace Scine
