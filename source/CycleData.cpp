#include "CycleData.h"
#include "Molecule.h"
#include "StdlibTypeAlgorithms.h"

#include <iostream>

namespace MoleculeManip {

/* Constructors */
CycleData::CycleData(const GraphType& sourceGraph) 
  : _baseGraphReference(sourceGraph)
{
  // Initialize
  _graphPtr = RDL_initNewGraph(boost::num_vertices(sourceGraph));

  // Copy vertices (without bond type information)
  auto iterPair = boost::edges(sourceGraph);
  for(auto iter = iterPair.first; iter != iterPair.second; iter++) {
    RDL_addUEdge(
      _graphPtr,
      boost::source(*iter, sourceGraph),
      boost::target(*iter, sourceGraph)
    );
  }

  // Calculate
  _dataPtr = RDL_calculate(_graphPtr);

  // Ensure it works
  assert(_dataPtr != nullptr);

}

/* Destructor */
CycleData::~CycleData() {
  RDL_deleteData(_dataPtr); // this also deletes the graph
}

/* Information */
unsigned CycleData::numCycleFamilies() const {
  return RDL_getNofURF(_dataPtr);
}

unsigned CycleData::numCycleFamilies(const AtomIndexType& index) const {
  return RDL_getNofURFContainingNode(_dataPtr, index);
}

//! Returns the number of relevant cycles (RCs)
unsigned CycleData::numRelevantCycles() const {
  return RDL_getNofRC(_dataPtr);
}

unsigned CycleData::numRelevantCycles(const AtomIndexType& index) const {
  return RDL_getNofRCFContainingNode(_dataPtr, index);
}

RDL_data* CycleData::getDataPtr() {
  return _dataPtr;
}

CycleIterator::CycleIterator(
  const CycleData& cycleData,
  const boost::optional<unsigned>& maxCycleSizeOption,
  const boost::optional<AtomIndexType>& containingIndexOption
) : _cycleDataRef(cycleData),
    _maxCycleSizeOption(maxCycleSizeOption),
    _containingIndexOption(containingIndexOption)
{
  // Initialize the cycle iterator
  _cycleIteratorPtr = RDL_getRCyclesIterator(cycleData._dataPtr);

  // Skip to first permissible ring per cycle size restrictions
  /* NOTE: While-do construction here because if the current cycle is
   * permisisble, we do not want to advance the cycle iterator
   */
  while(
    !atEnd()
    && !_currentCyclePermissible() 
  ) {
    RDL_cycleIteratorNext(_cycleIteratorPtr);
  }
}

bool CycleIterator::atEnd() const {
  return RDL_cycleIteratorAtEnd(_cycleIteratorPtr) != 0;
}

void CycleIterator::advance() {
  assert(!atEnd());
  /* Do-while construction to ensure that the iterator gets advance and only
   * after advancing it is checked whether the new cycle is permissible. 
   */
  do {
    RDL_cycleIteratorNext(_cycleIteratorPtr);
  } while(
    !atEnd()
    && !_currentCyclePermissible() 
  );
}

bool CycleIterator::_currentCyclePermissible() const {
  // Are there even constraints?
  if(_maxCycleSizeOption == boost::none && _containingIndexOption == boost::none) {
    return true;
  }

  RDL_cycle* cyclePtr = RDL_cycleIteratorGetCycle(_cycleIteratorPtr);
  unsigned cycleSize = cyclePtr -> weight;

  // Test against cycle size
  if(
    _maxCycleSizeOption 
    && cycleSize > _maxCycleSizeOption.value()
  ) {
    RDL_deleteCycle(cyclePtr);
    return false;
  }

  // Test if it contains the index 
  if(_containingIndexOption) {
    bool containsIndex = false;

    for(unsigned i = 0; i < cyclePtr -> weight; i++) {
      if(
        _containingIndexOption.value() == cyclePtr -> edges[i][0]
        || _containingIndexOption.value() == cyclePtr -> edges[i][1]
      ) {
        containsIndex = true;
        break;
      }
    }

    if(!containsIndex) {
      RDL_deleteCycle(cyclePtr);
      return false;
    }
  }

  RDL_deleteCycle(cyclePtr);
  return true;
}

CycleIterator::EdgeSet CycleIterator::getCurrentCycle() const {
  EdgeSet cycleEdges;

  assert(!atEnd());

  if(!atEnd()) {
    RDL_cycle* cyclePtr = RDL_cycleIteratorGetCycle(_cycleIteratorPtr);

    // Collect edge_descriptors for every edge of this cycle
    for(unsigned i = 0; i < cyclePtr -> weight; i++) {
      auto edgeIndexFoundPair = boost::edge(
        cyclePtr -> edges[i][0],
        cyclePtr -> edges[i][1],
        _cycleDataRef._baseGraphReference
      );

      /* These edges HAVE to be found, there should be a 1:1 mapping between
       * vertex indices
       */
      assert(edgeIndexFoundPair.second);

      cycleEdges.insert(edgeIndexFoundPair.first);
    }

    RDL_deleteCycle(cyclePtr);
  }

  return cycleEdges;
}

unsigned CycleIterator::cycleSize() const {
  assert(!atEnd());
  RDL_cycle* cyclePtr = RDL_cycleIteratorGetCycle(_cycleIteratorPtr);
  unsigned cycleSize = cyclePtr -> weight;
  RDL_deleteCycle(cyclePtr);
  return cycleSize;
}

CycleIterator::~CycleIterator() {
  RDL_deleteCycleIterator(_cycleIteratorPtr);
}

CycleIterator CycleData::getCyclesIterator() const {
  return CycleIterator(*this);
}

CycleIterator CycleData::getCyclesIteratorSizeLE(
  const unsigned& maxCycleSize
) const {
  return CycleIterator(*this, maxCycleSize, boost::none);
}

CycleIterator CycleData::getCyclesIteratorContaining(
  const AtomIndexType& containingIndex
) const {
  return CycleIterator(*this, boost::none, containingIndex);
}

std::map<AtomIndexType, unsigned> makeSmallestCycleMap(
  const CycleData& cycleData,
  const GraphType& graph
) {
  std::map<AtomIndexType, unsigned> smallestCycle;

  for(
    auto cycleIter = cycleData.getCyclesIterator();
    !cycleIter.atEnd();
    cycleIter.advance()
  ) {
    const auto cycleEdges = cycleIter.getCurrentCycle();
    const unsigned cycleSize = cycleEdges.size();

    for(const auto& edge : cycleEdges) {
      std::array<AtomIndexType, 2> indices {
        boost::source(edge, graph),
        boost::target(edge, graph)
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

std::vector<AtomIndexType> makeRingIndexSequence(
  const CycleIterator::EdgeSet& edgeSet,
  const GraphType& graph
) {
  // copy the edges so we can modify
  auto edgeDescriptors = edgeSet;

  auto firstEdgeIter = edgeDescriptors.begin();
  // Initialize with last edge descriptor's indices
  std::vector<AtomIndexType> indexSequence {
    boost::source(*firstEdgeIter, graph),
    boost::target(*firstEdgeIter, graph)
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
      if(boost::source(edge, graph) == indexSequence.back()) {
        indexSequence.emplace_back(
          boost::target(edge, graph)
        );
        edgeDescriptors.erase(edgeIter);
        break;
      }
      
      if(boost::target(edge, graph) == indexSequence.back()) {
        indexSequence.emplace_back(
          boost::source(edge, graph)
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


std::vector<AtomIndexType> centralizeRingIndexSequence(
  std::vector<AtomIndexType> ringIndexSequence,
  const AtomIndexType& center
) {
  assert(ringIndexSequence.front() == ringIndexSequence.back());

  ringIndexSequence.erase(ringIndexSequence.end() - 1);

  auto findIter = std::find(
    ringIndexSequence.begin(),
    ringIndexSequence.end(),
    center
  );

  assert(findIter != ringIndexSequence.end());

  std::rotate(
    ringIndexSequence.begin(),
    findIter,
    ringIndexSequence.end()
  );

  ringIndexSequence.push_back(
    center
  );

  return ringIndexSequence;
}


unsigned countPlanarityEnforcingBonds(
  const CycleIterator::EdgeSet& edgeSet,
  const GraphType& graph
) {
  unsigned count = 0;

  for(const auto& edge: edgeSet) {
    const auto& edgeType = graph[edge].bondType;
    if(
      edgeType == BondType::Double
      || edgeType == BondType::Aromatic
    ) {
      count += 1;
    }
  }

  return count;
}

} // namespace MoleculeManip
