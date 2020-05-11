/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Descriptors.h"

#include "molassembler/Graph.h"
#include "molassembler/Cycles.h"
#include "molassembler/Molecule.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/BondStereopermutator.h"

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

} // namespace molassmembler
} // namespace Scine
