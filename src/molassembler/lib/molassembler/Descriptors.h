#ifndef INCLUDE_MOLASSEMBLER_DESCRIPTORS_H
#define INCLUDE_MOLASSEMBLER_DESCRIPTORS_H

#include "Molecule.h"

namespace molassembler {

/*! Calculates number of freely rotatable bonds in a molecule
 *
 * Criteria:
 * - Must be a single bond
 * - Neither atom connected by the bond may be terminal
 * - May not be member of a small cycle. The threshold of this criterion
 *   can be set by the user. An edge that is a member of a cycle smaller than
 *   the threshold is not considered rotatable. The default is 5.
 *
 */
unsigned numRotatableBonds(
  const Molecule& mol,
  const unsigned cycleThreshold = 5
) {
  CycleData cycleData = mol.getCycleData();

  std::map<GraphType::edge_descriptor, unsigned> smallestCycle;

  for(
    auto cycleIter = cycleData.getCyclesIterator();
    !cycleIter.atEnd();
    cycleIter.advance()
  ) {
    const auto cycleEdges = cycleIter.getCurrentCycle();
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

} // namespace molassembler

#endif
