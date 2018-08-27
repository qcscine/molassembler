#include "molassembler/Descriptors.h"

#include "boost/range/iterator_range_core.hpp"
#include "temple/Containers.h"

#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/Cycles.h"
#include "molassembler/Molecule.h"

namespace molassembler {

unsigned numRotatableBonds(
  const Molecule& mol,
  const unsigned cycleThreshold
) {
  Cycles cycleData = mol.cycles();

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
          return cycleSize < currentMinCycleSize;
        }
      );
    }
  }

  unsigned count = 0;
  for(const auto& edge : boost::make_iterator_range(mol.graph().bonds())) {
    if(
      mol.graph().bondType(edge) == BondType::Single
      && (smallestCycle.count(edge) == 0 || smallestCycle.at(edge) > cycleThreshold)
    ) {
      // Neither of the atoms connected by this edge may be terminal
      if(
        mol.getNumAdjacencies(edge.first) > 1
        && mol.getNumAdjacencies(edge.second) > 1
      ) {
        ++count;
      }
    }
  }

  return count;
}

} // namespace molassmembler
