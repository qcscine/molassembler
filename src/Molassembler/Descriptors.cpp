/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Descriptors.h"

#include "Molassembler/Graph.h"
#include "Molassembler/Cycles.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/Graph/GraphAlgorithms.h"

#include "boost/dynamic_bitset.hpp"

namespace Scine {
namespace Molassembler {

unsigned numRotatableBonds(const Molecule& mol) {
  Cycles cycleData = mol.graph().cycles();

  std::unordered_map<BondIndex, unsigned, boost::hash<BondIndex>> smallestCycle;

  for(const auto& cycleEdges : cycleData) {
    const unsigned cycleSize = cycleEdges.size();

    for(const BondIndex& edge : cycleEdges) {
      auto findIter = smallestCycle.find(edge);

      if(findIter != std::end(smallestCycle)) {
        if(cycleSize < findIter->second) {
          findIter->second = cycleSize;
        }
      } else {
        smallestCycle.emplace(edge, cycleSize);
      }
    }
  }

  double count = 0;
  for(const auto& edge : mol.graph().bonds()) {
    // If the bond is not Single, it cannot be a rotatable bond
    if(mol.graph().bondType(edge) != BondType::Single) {
      continue;
    }

    // If either atom on the bond is terminal, it cannot be rotatable
    if(
      mol.graph().degree(edge.first) == 1
      || mol.graph().degree(edge.second) == 1
    ) {
      continue;
    }

    // If there is an assigned stereopermutator on the edge, it cannot be rotatable
    auto bondStereopermutatorOption = mol.stereopermutators().option(edge);
    if(
      bondStereopermutatorOption
      && bondStereopermutatorOption->assigned() != boost::none
    ) {
      continue;
    }

    auto findIter = smallestCycle.find(edge);

    // Is the edge part of a cycle?
    if(findIter == std::end(smallestCycle)) {
      // If not, the edge counts as a whole rotatable bond
      count += 1.0;
    } else {
      /* Otherwise, the contribution to the number of rotatable bonds is
       * calculated as
       *
       *   max(0.0, (cycle size - 3) / cycle size)
       *
       */
      unsigned cycleSize = findIter->second;
      count += std::max(0.0, (cycleSize - 3.0) / cycleSize);
    }
  }

  return static_cast<unsigned>(
    std::round(count)
  );
}

std::vector<AtomIndex> nonRankingEquivalentAtoms(const Molecule& mol) {
  std::unordered_set<AtomIndex> nonRankingEquivalents;
  for(const AtomIndex i : mol.graph().atoms()) {
    nonRankingEquivalents.insert(i);
  }

  /* Since atom stereopermutators are placed in an unordered_map, but the
   * following algorithm is sequence-dependent on them (a different ordering
   * yields equivalent sets of non-ranking-equivalent atoms, but different
   * ones), we sequence the permutator ordering here.
   */
  std::vector<AtomIndex> permutatorPlacements;
  permutatorPlacements.reserve(mol.stereopermutators().A());
  for(
    const AtomStereopermutator& permutator:
    mol.stereopermutators().atomStereopermutators()
  ) {
    permutatorPlacements.push_back(permutator.placement());
  }
  std::sort(std::begin(permutatorPlacements), std::end(permutatorPlacements));

  for(const AtomIndex placement : permutatorPlacements) {
    const auto& permutator = mol.stereopermutators().at(placement);

    // Do not try to get any information from vertices already marked duplicate
    if(nonRankingEquivalents.count(placement) == 0) {
      continue;
    }

    // Skip equivalently ranked substituents of stereogenic permutators
    if(permutator.numStereopermutations() > 1) {
      continue;
    }

    for(
      const auto& equallyRankedSubstituentGroup:
      permutator.getRanking().substituentRanking
    ) {
      // Skip single-vertex groups, those are uninteresting, really
      if(equallyRankedSubstituentGroup.size() == 1) {
        continue;
      }

      constexpr AtomIndex split = std::numeric_limits<AtomIndex>::max();
      const auto descendants = GraphAlgorithms::bfsUniqueDescendants(
        placement,
        equallyRankedSubstituentGroup,
        mol.graph().inner()
      );

      const auto groupBegin = std::begin(equallyRankedSubstituentGroup);
      const auto groupEnd = std::end(equallyRankedSubstituentGroup);

      for(const AtomIndex i : mol.graph().atoms()) {
        const AtomIndex descendant = descendants.at(i);
        // The split group of vertices is set as equivalent
        if(descendant == split) {
          nonRankingEquivalents.erase(i);
          continue;
        }

        /* If the vertex is a descendant of the equally ranked vertices
         * with the exception of the first equally ranked vertex (so as to keep
         * one of X equivalent branches), mark the vertex ranking equivalent
         */
        const auto findIter = std::find(groupBegin, groupEnd, descendant);
        if(findIter != groupBegin && findIter != groupEnd) {
          nonRankingEquivalents.erase(i);
        }
      }
    }
  }

  return std::vector<AtomIndex>(
    std::begin(nonRankingEquivalents),
    std::end(nonRankingEquivalents)
  );
}

} // namespace Molassembler
} // namespace Scine
