#include "DistanceGeometry/BFSConstraintCollector.h"

#include "CommonTrig.h"
#include "CNStereocenter.h"
#include "EZStereocenter.h"
#include "StdlibTypeAlgorithms.h"
#include "Log.h"

#include "cyclic_polygons/Minimal.h"

namespace MoleculeManip {

namespace DistanceGeometry {

/* Static public constexpr members declaration required for address access 
 * during runtime. Would exist compile-time only without
 */
constexpr double BFSConstraintCollector::bondRelativeVariance;
constexpr double BFSConstraintCollector::angleAbsoluteVariance;

BFSConstraintCollector::BFSConstraintCollector(
  const AdjacencyList& adjacencies,
  const StereocenterList& stereocenterList,
  MoleculeSpatialModel& spatialModel,
  const DistanceMethod& distanceMethod
) : _adjacencies(adjacencies),
    _distanceMethod(distanceMethod),
    _spatialModel(spatialModel),
    _cycleData(adjacencies.access()), // Construct cycleData from the graph
    _smallestCycleMap(_makeSmallestCycleMap()) // Construct the cycle size map
{
  // Check constraints on static constants
  static_assert(
    0 < bondRelativeVariance && bondRelativeVariance < 1,
    "BFSConstraintCollector static constant bond relative variance must fulfill"
    "0 < x << 1"
  );
  /* Unfortunately, this cannot be static_assert (although it is compile-time
   * constant due to C++14 limitations, see symmetry_information library for
   * possible C++17 improvements.
   */
  assert(
    0 < angleAbsoluteVariance && angleAbsoluteVariance < Symmetry::smallestAngle
    && "BFSConstraintCollector static constant angle absolute variance must fulfill"
    "0 < x << (smallest angle any local symmetry returns)"
  );

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
        // Instantiate without regard for number of assignments
        auto newStereocenterPtr = std::make_shared<Stereocenters::EZStereocenter>(
          source,
          adjacencies.rankPriority(target, {source}), // exclude shared edge
          target,
          adjacencies.rankPriority(source, {target})
        );

        // Map source *and* target to the same stereocenterPtr
        _stereocenterMap[source] = newStereocenterPtr;
        _stereocenterMap[target] = newStereocenterPtr;
      }
    }
  }

  // TODO these need tolerance adaptation in dependence of cycles too!
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
        _spatialModel.setDihedralBoundsIfEmpty(
          dihedralLimit.indices,
          dihedralLimit.lower,
          dihedralLimit.upper
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
      _stereocenterMap[i] = std::make_shared<Stereocenters::CNStereocenter>(
        adjacencies.determineLocalGeometry(i),
        i,
        adjacencies.rankPriority(i)
      );

      // At this point, new stereocenters should have only one assignment
      assert(_stereocenterMap[i] -> numAssignments() == 1);

      // Default assign it to zero
      _stereocenterMap[i] -> assign(0);
    }
  }

  /* For all cycles, set internal angles
   */
  for(
    auto cycleIter = _cycleData.getCyclesIteratorSizeLE(5);
    !cycleIter.atEnd();
    cycleIter.advance()
  ) {
    unsigned cycleSize = cycleIter.cycleSize();
    /* Gather sequence of atoms in cycle by progressively converting edge
     * descriptors into vertex indices
     */
    auto edgeDescriptors = cycleIter.getCurrentCycle();
    auto indexSequence = _makeRingIndexSequence(edgeDescriptors);

    /* There are a variety of cases here which all need to be treated
     * differently:
     * - Cycle of size 3: Always flat, set angles from cyclic polygons library
     * - Cycle of size 4: If flat (maybe already with a single double bond), use
     *   library to get precise internal angles and set them. Otherwise not all
     *   vertices are coplanar and all angles and dihedrals get a general
     *   tolerance increase
     * - Cycle of size 5: If aromatic, flat and we get precise internal angles
     *   from library. Otherwise, slightly increased angle and dihedral
     *   tolerances
     *
     * Assumptions
     * - Need to find out whether all vertices in cycles of size four are
     *   coplanar if one double bond is present. That seems to be the case for
     *   the few strained molecules collected in tests/strained_... Some more
     *   of these size four cycles are coplanar, especially when heteroatoms
     *   are present, but I don't immediately see the distinguishing criterion.
     */
    if( // flat cases
      cycleSize == 3
      || (
        cycleSize == 4
        && _countPlanarityEnforcingBonds(edgeDescriptors) >= 1
      )
      /* TODO missing cases: 
       * - Need aromaticity checking routine for cycle size 5
       */
    ) {
      /* First, we fetch the angles that maximize the cycle area using the
       * cyclic polygon library.
       */
      auto cycleInternalAngles = CyclicPolygons::internalAngles(
        // Map sequential index pairs to their purported bond length
        TemplateMagic::pairwiseMap( 
          indexSequence,
          [&](const AtomIndexType& i, const AtomIndexType& j) -> double {
            auto bondTypeOption = adjacencies.getBondType(i, j);

            // These vertices really ought to be bonded
            assert(bondTypeOption);

            return Bond::calculateBondDistance(
              adjacencies.getElementType(i),
              adjacencies.getElementType(j),
              bondTypeOption.value()
            );
          }
        )
      );

      /* The first angle returned is the angle between edges one and two, which
       * are indices [0, 1] and [1, 2] respectively. The last angle is between
       * edges [n - 1, 0], [0, 1].
       *
       * So, e.g. C4, C7, N9, O3::
       *
       *   index sequence 4, 7, 9, 3, 4
       *   distances  C-C, C-N, N-O, O-C
       *   angles  C-C-N, C-N-O, N-O-C, O-C-C
       *                                ^- Last overlaps past index sequence
       *
       * Now we set all internal angles with low tolerance
       */
      for(
        unsigned angleCentralIndex = 1;
        angleCentralIndex < indexSequence.size() - 1;
        ++angleCentralIndex
      ) {
        _spatialModel.setAngleBoundsIfEmpty(
          {
            indexSequence.at(angleCentralIndex - 1),
            indexSequence.at(angleCentralIndex),
            indexSequence.at(angleCentralIndex + 1),
          },
          cycleInternalAngles.at(angleCentralIndex - 1),
          BFSConstraintCollector::angleAbsoluteVariance
        );
        
        _angleSequences.insert(std::vector<AtomIndexType>{
          indexSequence.at(angleCentralIndex - 1),
          indexSequence.at(angleCentralIndex),
          indexSequence.at(angleCentralIndex + 1)
        });
      }
    } else { // Non-flat cases, tolerance increases on angles (and dihedrals)
    }
  }
  
}

std::vector<AtomIndexType> BFSConstraintCollector::_makeRingIndexSequence(
  const std::set<GraphType::edge_descriptor>& edgeSet
) const {
  // copy the edges so we can modify
  auto edgeDescriptors = edgeSet;

  auto firstEdgeIter = edgeDescriptors.begin();
  // Initialize with last edge descriptor's indices
  std::vector<AtomIndexType> indexSequence {
    boost::source(*firstEdgeIter, _adjacencies.access()),
    boost::target(*firstEdgeIter, _adjacencies.access())
  };
  edgeDescriptors.erase(firstEdgeIter);
  // firstEdgeIter is now invalid!

  while(!edgeDescriptors.empty()) {
    for(
      auto edgeIter = edgeDescriptors.begin();
      edgeIter != edgeDescriptors.end();
      ++edgeIter
    ) {
      auto& edge = *edgeIter;
      if(boost::source(edge, _adjacencies.access()) == indexSequence.back()) {
        indexSequence.emplace_back(
          boost::target(edge, _adjacencies.access())
        );
        edgeDescriptors.erase(edgeIter);
        break;
      }
      
      if(boost::target(edge, _adjacencies.access()) == indexSequence.back()) {
        indexSequence.emplace_back(
          boost::source(edge, _adjacencies.access())
        );
        edgeDescriptors.erase(edgeIter);
        break;
      }
    }
  }
  /* Now indexSequence should contain the entire sequence, but the first
   * vertex index occurs twice, once at the front and once at the back. 
   */

  return indexSequence;
}

unsigned BFSConstraintCollector::_countPlanarityEnforcingBonds(
  const std::set<GraphType::edge_descriptor>& edgeSet
) const {
  unsigned count = 0;

  for(const auto& edge: edgeSet) {
    const auto& edgeType = _adjacencies.access()[edge].bondType;
    if(
      edgeType == BondType::Double
      || edgeType == BondType::Aromatic
    ) {
      count += 1;
    }
  }

  return count;
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

      _spatialModel.setBondBounds(
        {i, j},
        bondDistance,
        BFSConstraintCollector::bondRelativeVariance
      );
    }
  } else { // Uniform
    constexpr double lower = (1 - BFSConstraintCollector::bondRelativeVariance) * 1.5;
    constexpr double upper = (1 + BFSConstraintCollector::bondRelativeVariance) * 1.5;

    for(const auto& edge: _adjacencies.iterateEdges()) {
      unsigned i = boost::source(edge, _adjacencies.access());
      unsigned j = boost::target(edge, _adjacencies.access());

      _spatialModel.setBondBounds(
        {i, j},
        lower,
        upper
      );
    }
  }
}

void BFSConstraintCollector::set13Bounds(const std::vector<AtomIndexType>& chain) {
  assert(chain.size() == 3);

  _spatialModel.setAngleBoundsIfEmpty(
    {chain.front(), chain.at(1), chain.back()},
    _getAngle(
      chain.front(),
      chain.at(1),
      chain.back()
    ),
    BFSConstraintCollector::angleAbsoluteVariance
  );
}

void BFSConstraintCollector::set14Bounds(const std::vector<AtomIndexType>& chain) {
  /* TODO This is largely nonsense
   * - Now, it would be better to just linearly extract all angle information
   *   from the Stereocenters and dihedral information from EZStereocenters. 
   *   Apart from that, all other possible dihedrals are obviously restricted
   *   to [0, π], but there's no reason to BFSVisit everything for that
   */
  _spatialModel.setDihedralBoundsIfEmpty(
    {
      chain.front(),
      chain.at(1),
      chain.at(2),
      chain.back()
    },
    0,
    M_PI
  );
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

BFSConstraintCollector::SmallestCycleMapType BFSConstraintCollector::_makeSmallestCycleMap() const {
  SmallestCycleMapType smallestCycle;

  for(
    auto cycleIter = _cycleData.getCyclesIterator();
    !cycleIter.atEnd();
    cycleIter.advance()
  ) {
    const auto cycleEdges = cycleIter.getCurrentCycle();
    const unsigned cycleSize = cycleEdges.size();

    for(const auto& edge : cycleEdges) {
      std::array<AtomIndexType, 2> indices {
        boost::source(edge, _adjacencies.access()),
        boost::target(edge, _adjacencies.access())
      };

      for(const auto& index: indices) {
        StdlibTypeAlgorithms::addOrUpdateMapIf(
          smallestCycle,
          index, // key_type to check
          cycleSize, // mapped_value to place if key does not exist or if ...
          [&cycleSize](const unsigned& currentMinCycleSize) -> bool {
            return cycleSize < currentMinCycleSize;
          }
        );
      }
    }
  }

  return smallestCycle;
}

} // namespace DistanceGeometry

} // namespace MoleculeManip
