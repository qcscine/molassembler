#include "molassembler/Descriptors.h"

#include "boost/range/iterator_range_core.hpp"

#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/Cycles.h"
#include "molassembler/Molecule.h"
#include "molassembler/StereocenterList.h"

namespace molassembler {

unsigned numRotatableBonds(const Molecule& mol) {
  Cycles cycleData = mol.graph().cycles();

  std::map<BondIndex, unsigned> smallestCycle;

  for(const auto cyclePtr : cycleData) {
    const auto cycleEdges = Cycles::edges(cyclePtr);
    const unsigned cycleSize = cycleEdges.size();

    for(const auto& edge : cycleEdges) {
      StdlibTypeAlgorithms::addOrUpdateMapIf(
        smallestCycle,
        edge,
        cycleSize,
        [&cycleSize](const unsigned currentMinCycleSize) -> bool {
          return cycleSize <= currentMinCycleSize;
        }
      );
    }
  }

  double count = 0;
  for(const auto& edge : boost::make_iterator_range(mol.graph().bonds())) {
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

    // If there is an assigned stereocenter on the edge, it cannot be rotatable
    auto bondStereocenterOption = mol.stereocenters().option(edge);
    if(
      bondStereocenterOption
      && bondStereocenterOption->assigned() != boost::none
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

  return std::round(count);
}

} // namespace molassmembler
