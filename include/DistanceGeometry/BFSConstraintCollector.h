#ifndef INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H
#define INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H

#include "GraphDistanceMatrix.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"

namespace MoleculeManip {

namespace DistanceGeometry {

class BFSConstraintCollector {
private:
  /* Private Members */
  // input
  const AdjacencyList& _adjacencies;

  // output (via side effect)
  DistanceBoundsMatrix& _distanceBounds;

  // mutable state
  AtomIndexType _initial;
  std::vector<
    std::vector<AtomIndexType>
  > _chains;

public:
  /* Ctor */
  BFSConstraintCollector(
    const AdjacencyList& adjacencies,
    DistanceBoundsMatrix& distanceBounds,
    const AtomIndexType& initial
  ) : _adjacencies(adjacencies),
      _distanceBounds(distanceBounds),
      _initial(initial)
  {}

  /* Function operator for BFSVisit, is impure */
  bool operator() (const AtomIndexType& atomIndex) {
    // continue here with constraint generation
    if(atomIndex == _initial) {
      // create an initial chain
      _chains.emplace_back(
        std::vector<AtomIndexType>({atomIndex})
      );
    } else {
      // add atomIndex to all chains that have your adjacencies at the back
    }

    // continue BFS
    return true;
  }

  /* To reuse a single ConstraintCollector, use this to reset internal state
   * for another starting atom index.
   */
  void resetToInitial(const AtomIndexType& initial) {
    _chains.clear();
    _initial = initial;
  }
};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
