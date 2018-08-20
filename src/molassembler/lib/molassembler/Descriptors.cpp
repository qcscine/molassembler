#include "molassembler/Descriptors.h"

#include "temple/Containers.h"

#include "molassembler/detail/StdlibTypeAlgorithms.h"
#include "molassembler/Cycles.h"
#include "molassembler/Molecule.h"

namespace molassembler {

unsigned numRotatableBonds(
  const Molecule& mol,
  const unsigned cycleThreshold
) {
  Cycles cycleData = mol.getCycleData();

  std::map<GraphType::edge_descriptor, unsigned> smallestCycle;

  for(const auto cyclePtr : cycleData) {
    const auto cycleEdges = Cycles::edges(cyclePtr, mol.getGraph());
    const unsigned cycleSize = cycleEdges.size();

    for(const auto& edge : cycleEdges) {
      StdlibTypeAlgorithms::addOrUpdateMapIf(
        smallestCycle,
        edge,
        cycleSize,
        [&cycleSize](const unsigned& currentMinCycleSize) -> bool {
          return cycleSize < currentMinCycleSize;
        }
      );
    }
  }

  unsigned count = 0;
  for(const auto& edge : mol.iterateEdges()) {
    if(
      mol.getBondType(edge) == BondType::Single
      && (smallestCycle.count(edge) == 0 || smallestCycle.at(edge) > cycleThreshold)
    ) {
      // Neither of the atoms connected by this edge may be terminal
      if(
        temple::all_of(
          mol.vertices(edge),
          [&mol](const auto& vertexIndex) -> bool {
            return mol.getNumAdjacencies(vertexIndex) > 1;
          }
        )
      ) {
        ++count;
      }
    }
  }

  return count;
}

} // namespace molassmembler
