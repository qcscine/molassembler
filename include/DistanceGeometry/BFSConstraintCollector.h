#ifndef INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H
#define INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H

#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "AdjacencyList.h"
#include "Tree.h"

#include "CommonTrig.h"
#include "CNStereocenter.h"
#include "EZStereocenter.h"

/* TODO
 * - how to get angles at atoms when no stereocenters exist?
 *   -> If there is a Stereocenter centered on the intermediate atom
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

  // Pre-set constants
  const double oneTwoVariance = 0.05;

  std::map<
    AtomIndexType,
    std::shared_ptr<Stereocenters::Stereocenter>
  > _stereocenterMap;

public:
  using NodeType = Tree::Node<AtomIndexType>;

  /* Ctor */
  /*! 
   * Assumes the Distance Bounds is entirely uninitialized (lower diagonal is 0,
   * upper diagonal is 100).
   */
  BFSConstraintCollector(
    const AdjacencyList& adjacencies,
    const StereocenterList& stereocenterList,
    DistanceBoundsMatrix& distanceBounds
  ) : _adjacencies(adjacencies),
      _distanceBounds(distanceBounds)
  {
    // Pre-fill the Distance bounds with proper 1-2 distances
    for(const auto& explicitEdgePair : adjacencies.getEdges()) {
      auto bondDistance = Bond::calculateBondDistance(
        adjacencies.getElementType(
          explicitEdgePair.first.first
        ),
        adjacencies.getElementType(
          explicitEdgePair.first.second
        ),
        explicitEdgePair.second
      );

      _distanceBounds.setLowerBound(
        explicitEdgePair.first.first,
        explicitEdgePair.first.second,
        bondDistance - oneTwoVariance
      );

      _distanceBounds.setUpperBound(
        explicitEdgePair.first.first,
        explicitEdgePair.first.second,
        bondDistance + oneTwoVariance
      );
    }

    // Set lower bounds to sum of vdw radii
    for(unsigned i = 0; i < adjacencies.numAtoms(); i++) {
      for(unsigned j = i + 1; j < adjacencies.numAtoms(); j++) {
        /* setting the bounds will fail for bonded pairs as those have strict
         * bounds already and the fairly high sum of vdw would lead to
         * inconsistencies
         */
        _distanceBounds.setLowerBound(i, j,
          AtomInfo::vdwRadius(
            adjacencies.getElementType(i)
          ) + AtomInfo::vdwRadius(
            adjacencies.getElementType(j)
          ) 
        );
      }
    }

    // Populate the stereocenterMap with copies of the molecule's stereocenters
    // Start with the existing stereocenters with multiple assignments
    for(const auto& stereocenterPtr : stereocenterList) {
      for(const auto& involvedAtom : stereocenterPtr -> involvedAtoms()) {
        if(_stereocenterMap.count(involvedAtom) == 0) {
          if(stereocenterPtr -> type() == Stereocenters::Type::CNStereocenter) {
            // Downcast the shared ptr
            auto CNSPtr = std::dynamic_pointer_cast<
              Stereocenters::CNStereocenter
            >(stereocenterPtr);

            // Explicit new shared ptr using derived copy-constructor
            _stereocenterMap[involvedAtom] = std::make_shared<
              Stereocenters::CNStereocenter
            >(*CNSPtr);
          } else {
            auto EZSPtr = std::dynamic_pointer_cast<
              Stereocenters::EZStereocenter
            >(stereocenterPtr);

            _stereocenterMap[involvedAtom] = std::make_shared<
              Stereocenters::EZStereocenter
            >(*EZSPtr);
          }

          /* Now we have a stereocenter, but it might be unassigned, in which
           * case angle() will fail! Need to assign unassigned ones
           * at random consistent with the unique assignments' relative
           * occurrences (TODO)
           */
          if(!_stereocenterMap[involvedAtom] -> assigned()) {
            _stereocenterMap[involvedAtom] -> assign(0);
          }
        }
      }
    }

    /* For every missing non-terminal atom, create a CNStereocenter in the
     * determined geometry
     */
    for(unsigned i = 0; i < adjacencies.numAtoms(); i++) {
      if(
        _stereocenterMap.count(i) == 0  // not already in the map
        && adjacencies.getNumAdjacencies(i) > 1 // non-terminal
      ) {
        auto rankResultPair = adjacencies.rankPriority(i);

        _stereocenterMap[i] = std::make_shared<Stereocenters::CNStereocenter>(
          adjacencies.determineLocalGeometry(i),
          i,
          rankResultPair.first,
          rankResultPair.second
        );

        assert(_stereocenterMap[i] -> assignments() == 1);

        // Default assign it to zero
        _stereocenterMap[i] -> assign(0);
      }
    }
  }

  double _toRadians(const double& inDegrees) {
    return M_PI * inDegrees / 180;
  }

  inline double _getAngle(
    AtomIndexType& i,
    AtomIndexType& j,
    AtomIndexType& k
  ) {
    return _toRadians(
      _stereocenterMap[j] -> angle(i, j, k)
    );
  }


  /* Function operator for Tree BFSVisit, is impure */
  bool operator() (
    const std::shared_ptr<NodeType>& nodePtr,
    const unsigned& depth
  ) {
    std::vector<AtomIndexType> chain;

    auto currentNode = nodePtr;
    while(!(currentNode -> parentWeakPtr).expired()) {
      chain.push_back(currentNode -> key);
      currentNode = currentNode -> parentWeakPtr.lock();
    }
    // Add root to the chain too, although it does not have a parent
    chain.push_back(currentNode -> key);

    if(depth == 2) { // angle
      double angle = _getAngle(
        chain.front(),
        chain.at(1),
        chain.back()
      );

      _distanceBounds.setLowerBound(
        chain.front(),
        chain.back(),
        CommonTrig::lawOfCosines(
          _distanceBounds.lowerBound(
            chain.front(),
            chain.at(1)
          ),
          _distanceBounds.lowerBound(
            chain.at(1),
            chain.back()
          ),
          angle
        )
      );

      _distanceBounds.setUpperBound(
        chain.front(),
        chain.back(),
        CommonTrig::lawOfCosines(
          _distanceBounds.upperBound(
            chain.front(),
            chain.at(1)
          ),
          _distanceBounds.upperBound(
            chain.at(1),
            chain.back()
          ),
          angle
        )
      );
    } else if(depth == 3) { // dihedral
      double abAngle = _getAngle(
        chain.front(),
        chain.at(1),
        chain.at(2)
      );
      double bcAngle = _getAngle(
        chain.at(1),
        chain.at(2),
        chain.back()
      ); 

      _distanceBounds.setLowerBound(
        chain.front(),
        chain.back(),
        CommonTrig::dihedralLength(
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
        )
      );

      _distanceBounds.setUpperBound(
        chain.front(),
        chain.back(),
        CommonTrig::dihedralLength(
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
        )
      );
    } 

    // continue BFS
    return true;
  }

  std::vector<
    Stereocenters::Stereocenter::ChiralityConstraintPrototype
  > getChiralityPrototypes() const {
    std::vector<
      Stereocenters::Stereocenter::ChiralityConstraintPrototype
    > prototypes;

    for(const auto& iterPair : _stereocenterMap) {
      auto chiralityPrototypes = iterPair.second -> chiralityConstraints();

      std::copy(
        chiralityPrototypes.begin(),
        chiralityPrototypes.end(),
        std::back_inserter(prototypes)
      );
    }

    return prototypes;
  }

// TODO: move this to some other place in the DG process,
// AFTER a distances matrix has been created, the prototypes can be altered to
// real constraints
//  std::vector<ChiralityConstraint> getChiralityConstraints(
//    const Eigen::MatrixXd& distancesMatrix
//  ) {
//  /*
//   * Return a list of chirality constraints.  The target volume of the
//   * chirality constraint created by the tetrahedron is calculated using
//   * internal coordinates (the Cayley-Menger determinant), always leading to V
//   * > 0, so depending on the current assignment, the sign of the result is
//   * switched. The formula used later in chirality constraint calculation for
//   * explicit coordinates is adjusted by V' = 6 V to avoid an unnecessary
//   * factor, so we do that here too:
//   *               
//   *    288 V²  = |...|               | substitute V' = 6 V
//   * -> 8 (V')² = |...|               
//   * ->      V' = sqrt(|...| / 8)
//   *
//   * where the Cayley-Menger determinant |...| is square symmetric:
//   *   
//   *          |   0    1    1    1    1  |
//   *          |        0  d12² d13² d14² |
//   *  |...| = |             0  d23² d24² |
//   *          |                  0  d34² |
//   *          |  ...                  0  |
//   *
//   */
//    std::vector<ChiralityConstraint> constraints;
//
//    for(const auto& iterPair : _stereocenterMap) {
//      auto chiralityPrototypes = iterPair.second -> chiralityConstraints();
//
//      for(const auto& prototype : chiralityPrototypes) {
//        Eigen::Matrix<double, 5, 5> cayleyMenger;
//        cayleyMenger.setZero();
//
//        for(unsigned i = 0; i < 4; i++) {
//          for(unsigned j = i + 1; j < 4; j++) {
//            cayleyMenger(j + 1, i + 1) = pow(
//              distancesMatrix(
//                std::max(
//                  prototype.first[i],
//                  prototype.first[j]
//                ),
//                std::min(
//                  prototype.first[i],
//                  prototype.first[j]
//                )
//              ),
//              2
//            );
//          }
//        }
//
//        // top row of cayleyMenger matrix
//        for(unsigned i = 1; i < 5; i++) {
//          cayleyMenger(0, i) = 1;
//        }
//
//        auto determinant = static_cast<
//          Eigen::Matrix<double, 5, 5>
//        >(
//          cayleyMenger.selfadjointView<Eigen::Upper>()
//        ).determinant();
//
//        auto chiralityTarget = sqrt(
//          determinant / 8.0
//        );
//
//        if(
//          prototype.second 
//          == Stereocenters::Stereocenter::ChiralityConstraintTarget::Negative
//        ) {
//          chiralityTarget *= -1;
//        }
//      }
//    }
//  }
};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
