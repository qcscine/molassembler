#include "DistanceGeometry/BFSConstraintCollector.h"

#include "CommonTrig.h"
#include "CNStereocenter.h"
#include "EZStereocenter.h"

#include "Log.h"


namespace MoleculeManip {

namespace DistanceGeometry {

BFSConstraintCollector::BFSConstraintCollector(
  const AdjacencyList& adjacencies,
  const StereocenterList& stereocenterList,
  DistanceBoundsMatrix& distanceBounds,
  const DistanceMethod& distanceMethod
) : _adjacencies(adjacencies),
    _distanceMethod(distanceMethod),
    _distanceBounds(distanceBounds)
{
  if(_distanceMethod == DistanceMethod::UFFLike) {
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
  } else { // Uniform
    for(const auto& explicitEdgePair : adjacencies.getEdges()) {
      _distanceBounds.setLowerBound(
        explicitEdgePair.first.first,
        explicitEdgePair.first.second,
        1.5 - oneTwoVariance
      );

      _distanceBounds.setUpperBound(
        explicitEdgePair.first.first,
        explicitEdgePair.first.second,
        1.5 + oneTwoVariance
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

      assert(_stereocenterMap[i] -> numAssignments() == 1);

      // Default assign it to zero
      _stereocenterMap[i] -> assign(0);
    }
  }
}

/* impure Function operator for Tree BFSVisit */
bool BFSConstraintCollector::operator() (
  const std::shared_ptr<NodeType>& nodePtr,
  const unsigned& depth
) {

#ifndef NDEBUG
  Log::log(Log::Particulars::BFSConstraintCollectorVisitCall) 
    << "BFSConstraintCollector: operator() call at depth " 
    << depth << " on node with key " << nodePtr -> key << std::endl;
#endif

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

#ifndef NDEBUG
    Log::log(Log::Particulars::BFSConstraintCollectorVisitCall) 
      << "BFSConstraintCollector: attempting to improve bounds on chain {" 
      << chain.front() << ", " << chain.at(1) << ", " << chain.back() << "}, "
      << "angle: " << angle << ", suggest bounds ["
      << CommonTrig::lawOfCosines(
        _distanceBounds.lowerBound(
          chain.front(),
          chain.at(1)
        ),
        _distanceBounds.lowerBound(
          chain.at(1),
          chain.back()
        ),
        angle
      ) << ", " << CommonTrig::lawOfCosines(
        _distanceBounds.upperBound(
          chain.front(),
          chain.at(1)
        ),
        _distanceBounds.upperBound(
          chain.at(1),
          chain.back()
        ),
        angle
      ) << "], currently [" 
      << _distanceBounds.lowerBound(chain.front(), chain.back())
      << ", " << _distanceBounds.upperBound(chain.front(), chain.back())
      << "]" << std::endl;
#endif
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

void BFSConstraintCollector::finalizeBoundsMatrix() {
  // Set lower bounds to sum of vdw radii
  for(unsigned i = 0; i < _adjacencies.numAtoms(); i++) {
    for(unsigned j = i + 1; j < _adjacencies.numAtoms(); j++) {
      /* setting the bounds will fail for bonded pairs as those have strict
       * bounds already and the fairly high sum of vdw would lead to
       * inconsistencies
       */
      if(_distanceBounds.lowerBound(i, j) == 0) {
        _distanceBounds.setLowerBound(i, j,
          AtomInfo::vdwRadius(
            _adjacencies.getElementType(i)
          ) + AtomInfo::vdwRadius(
            _adjacencies.getElementType(j)
          ) 
        );
      }
    }
  }

  _distanceBounds.smooth();

  assert(_distanceBounds.boundInconsistencies() == 0);
}

std::vector<
  Stereocenters::Stereocenter::ChiralityConstraintPrototype
> BFSConstraintCollector::getChiralityPrototypes() const {
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

} // namespace DistanceGeometry

} // namespace MoleculeManip
