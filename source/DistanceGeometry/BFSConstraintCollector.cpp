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
  set12Bounds();

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
         *
         * Make sure changes are only effected on stereocenters in our local
         * map of stereocenters, we are forbidden from changing the molecular
         * graph by altering assignments
         */
        if(!_stereocenterMap[involvedAtom] -> assigned()) {
          _stereocenterMap[involvedAtom] -> assign(0);
        }
      }
    }
  }

  // For every non-implicated double bond, create an EZStereocenter 
  for(const auto& edgeIndex: adjacencies.iterateEdges()) {
    if(adjacencies.access()[edgeIndex].bondType == BondType::Double) {
      auto source = boost::source(edgeIndex, adjacencies.access()),
           target = boost::target(edgeIndex, adjacencies.access());

      if(_stereocenterMap.count(source) == 0 && _stereocenterMap.count(target) == 0) {

        // Calculate Priorities for each's substituents
        auto sourceSubstituentsRanking = adjacencies.rankPriority(
          source,
          {target} // exclude edge sharing neighbor
        );

        auto targetSubstituentsRanking = adjacencies.rankPriority(
          target,
          {source} // exclude edge sharing neighbor
        );

        // Instantiate without regard for number of assignments
        auto newStereocenterPtr = std::make_shared<Stereocenters::EZStereocenter>(
          source,
          sourceSubstituentsRanking.first,
          sourceSubstituentsRanking.second,
          target,
          targetSubstituentsRanking.first,
          targetSubstituentsRanking.second
        );

        // Map source *and* target to the same stereocenterPtr
        _stereocenterMap[source] = newStereocenterPtr;
        _stereocenterMap[target] = newStereocenterPtr;
      }
    }
  }

  // For every EZStereocenter, get the dihedral information and set it
  for(const auto& mapIterPair : _stereocenterMap) {
    const auto& mappingIndex = mapIterPair.first;
    const auto& stereocenterPtr = mapIterPair.second;
    if(
      // Filter out CNStereocenters
      stereocenterPtr -> type() == Stereocenters::Type::EZStereocenter
      /* Since EZStereocenters are in the map twice (for each of their
       * involvedAtoms()), and we only want to process their dihedral
       * information once, check the mapping index against the first index in
       * the involvedAtoms set
       */
      && mappingIndex == *(
        stereocenterPtr -> involvedAtoms().begin()
      )
    ) {
      for(const auto& dihedralLimit : stereocenterPtr -> dihedralLimits()) {
        _distanceBounds.setLowerBound(
          dihedralLimit.indices.front(),
          dihedralLimit.indices.back(),
          CommonTrig::dihedralLength(
            _distanceBounds.lowerBound(
              dihedralLimit.indices.front(),
              dihedralLimit.indices.at(1)
            ),
            _distanceBounds.lowerBound(
              dihedralLimit.indices.at(1),
              dihedralLimit.indices.at(2)
            ),
            _distanceBounds.lowerBound(
              dihedralLimit.indices.at(2),
              dihedralLimit.indices.back()
            ),
            ConstexprMagic::Math::toRadians(120),
            ConstexprMagic::Math::toRadians(120),
            dihedralLimit.lower
          )
        );

        _distanceBounds.setUpperBound(
          dihedralLimit.indices.front(),
          dihedralLimit.indices.back(),
          CommonTrig::dihedralLength(
            _distanceBounds.upperBound(
              dihedralLimit.indices.front(),
              dihedralLimit.indices.at(1)
            ),
            _distanceBounds.upperBound(
              dihedralLimit.indices.at(1),
              dihedralLimit.indices.at(2)
            ),
            _distanceBounds.upperBound(
              dihedralLimit.indices.at(2),
              dihedralLimit.indices.back()
            ),
            ConstexprMagic::Math::toRadians(120),
            ConstexprMagic::Math::toRadians(120),
            dihedralLimit.upper
          )
        );

        /* Add the sequence considered to the dihedral sequences to prevent set
         * limits to be overwritten or altered further during tree traversals
         */
        _dihedralSequences.insert(dihedralLimit.indices);
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

      // At this point, new stereocenters should have only one assignment
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
    if(!_angleSequences.contains(chain)) {
      set13Bounds(chain);
      _angleSequences.insert(chain);
    }
  } else if(depth == 3) { // dihedral
    if(!_dihedralSequences.contains(chain)) {
      /* No need to check whether 2-3 is the same EZStereocenter, those sequences
       * are already in _dihedralSequences (see constructor)
       */
      set14Bounds(chain);
      _dihedralSequences.insert(chain);
    }
  }

  // continue BFS
  return true;
}

void BFSConstraintCollector::set12Bounds() {
  if(_distanceMethod == DistanceMethod::UFFLike) {
    for(const auto& edge: _adjacencies.iterateEdges()) {
      unsigned i = boost::source(edge, _adjacencies.access());
      unsigned j = boost::target(edge, _adjacencies.access());
      auto bondType = _adjacencies.access()[edge].bondType;

      auto bondDistance = Bond::calculateBondDistance(
        _adjacencies.getElementType(i),
        _adjacencies.getElementType(j),
        bondType
      );

      _distanceBounds.setLowerBound(
        i,
        j,
        bondDistance - oneTwoVariance
      );

      _distanceBounds.setUpperBound(
        i,
        j,
        bondDistance + oneTwoVariance
      );
    }
  } else { // Uniform
    for(const auto& edge: _adjacencies.iterateEdges()) {
      unsigned i = boost::source(edge, _adjacencies.access());
      unsigned j = boost::target(edge, _adjacencies.access());
      _distanceBounds.setLowerBound(i, j, 1.5 - oneTwoVariance);
      _distanceBounds.setUpperBound(i, j, 1.5 + oneTwoVariance);
    }
  }
}

void BFSConstraintCollector::set13Bounds(
  const std::vector<AtomIndexType>& chain
) {
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
}

void BFSConstraintCollector::set14Bounds(
  const std::vector<AtomIndexType>& chain
) {
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

void BFSConstraintCollector::finalizeBoundsMatrix() {
  assert(_distanceBounds.boundInconsistencies() == 0);

  // Set lower bounds to sum of vdw radii
  for(unsigned i = 0; i < _adjacencies.numAtoms(); i++) {
    for(unsigned j = i + 1; j < _adjacencies.numAtoms(); j++) {
      /* setting the bounds will fail for bonded pairs as those have strict
       * bounds already and the fairly high sum of vdw would lead to
       * inconsistencies
       */
      if(_distanceBounds.lowerBound(i, j) == 0) {
        _distanceBounds.setLowerBound(
          i,
          j,
          AtomInfo::vdwRadius(
            _adjacencies.getElementType(i)
          ) + AtomInfo::vdwRadius(
            _adjacencies.getElementType(j)
          ) 
        );
      }
    }
  }

  assert(_distanceBounds.boundInconsistencies() == 0);

  _distanceBounds.smooth();

  assert(_distanceBounds.boundInconsistencies() == 0);
}

std::vector<
  Stereocenters::ChiralityConstraintPrototype
> BFSConstraintCollector::getChiralityPrototypes() const {
  std::vector<
    Stereocenters::ChiralityConstraintPrototype
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
