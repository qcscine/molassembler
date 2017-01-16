#ifndef INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H
#define INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H

#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "AdjacencyList.h"
#include "Tree.h"
#include "TreeAlgorithms.h"
#include "CommonTrig.h"

/* TODO
 * - how to get angles at atoms when no stereocenters exist?
 * - tests
 */

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

public:
  using NodeType = Tree::Node<AtomIndexType>;

  /* Ctor */
  BFSConstraintCollector(
    const AdjacencyList& adjacencies,
    DistanceBoundsMatrix& distanceBounds,
    const AtomIndexType& initial
  ) : _adjacencies(adjacencies),
      _distanceBounds(distanceBounds),
      _initial(initial)
  {}

  /* Function operator for Tree BFSVisit, is impure */
  bool operator() (std::shared_ptr<NodeType>& nodePtr, unsigned depth) {
    std::vector<AtomIndexType> chain;
    for(
      std::shared_ptr<NodeType> currentNode = nodePtr;
      !(currentNode -> parentWeakPtr).expired();
      currentNode = currentNode -> parentWeakPtr.lock()
    ) {
      chain.push_back(currentNode -> key);
    }

    if(depth == 2) { // angle
      // TODO how to get angle?
      double angle = 0; 
      _distanceBounds.lowerBound(
        chain.front(),
        chain.back()
      ) = CommonTrig::lawOfCosines(
        _distanceBounds.lowerBound(
          chain.front(),
          chain.at(1)
        ),
        _distanceBounds.lowerBound(
          chain.at(1),
          chain.back()
        ),
        angle
      );
      _distanceBounds.upperBound(
        chain.front(),
        chain.back()
      ) = CommonTrig::lawOfCosines(
        _distanceBounds.upperBound(
          chain.front(),
          chain.at(1)
        ),
        _distanceBounds.upperBound(
          chain.at(1),
          chain.back()
        ),
        angle
      );
    } else if(depth == 3) { // dihedral

      // TODO get angles
      double abAngle = 0, bcAngle = 0; 

      _distanceBounds.lowerBound(
        chain.front(),
        chain.back()
      ) = CommonTrig::dihedralLength(
        _distanceBounds.lowerBound(
          chain.front(),
          chain.at(1)
        ),
        _distanceBounds.lowerBound(
          chain.at(1),
          chain.at(2)
        ),
        _distanceBounds.lowerBound(
          chain.at(2),
          chain.back()
        ),
        abAngle,
        bcAngle,
        0 // cis dihedral
      );

      _distanceBounds.upperBound(
        chain.front(),
        chain.back()
      ) = CommonTrig::dihedralLength(
        _distanceBounds.upperBound(
          chain.front(),
          chain.at(1)
        ),
        _distanceBounds.upperBound(
          chain.at(1),
          chain.at(2)
        ),
        _distanceBounds.upperBound(
          chain.at(2),
          chain.back()
        ),
        abAngle,
        bcAngle,
        M_PI // trans dihedral
      );
    } 

    // continue BFS
    return true;
  }

  /* To reuse a single ConstraintCollector, use this to reset internal state
   * for another starting atom index.
   */
  void resetToInitial(const AtomIndexType& initial) {
    _initial = initial;
  }
};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
