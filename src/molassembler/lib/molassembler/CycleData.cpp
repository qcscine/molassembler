#include "CycleData.h"
#include "StdlibTypeAlgorithms.h"

#include <iostream>

namespace molassembler {

Cycles::Cycles(const GraphType& sourceGraph) {
  _rdlPtr = std::make_shared<RDLDataPtrs>(sourceGraph);
}

unsigned Cycles::size(const RDL_cycle* const cyclePtr) {
  return cyclePtr -> weight;
}

Cycles::EdgeList Cycles::edgeVertices(const RDL_cycle* const cyclePtr) {
  EdgeList cycleEdges;

  // Collect the pair of atom indices for every edge of this cycle
  for(unsigned i = 0; i < size(cyclePtr); i++) {
    cycleEdges.emplace_back(
      std::array<AtomIndexType, 2> {
        cyclePtr->edges[i][0],
        cyclePtr->edges[i][1]
      }
    );
  }

  return cycleEdges;
}

temple::TinySet<GraphType::edge_descriptor> Cycles::edges(
  const RDL_cycle* const cyclePtr,
  const GraphType& graph
) {
  temple::TinySet<GraphType::edge_descriptor> cycleEdges;

  // Collect the pair of atom indices for every edge of this cycle
  for(unsigned i = 0; i < size(cyclePtr); ++i) {
    auto edgeFoundPair = boost::edge(
      cyclePtr->edges[i][0],
      cyclePtr->edges[i][1],
      graph
    );

    assert(edgeFoundPair.second);

    cycleEdges.insert(edgeFoundPair.first);
  }

  return cycleEdges;
}

temple::TinySet<GraphType::edge_descriptor> Cycles::edges(
  const EdgeList& edges,
  const GraphType& graph
) {
  temple::TinySet<GraphType::edge_descriptor> cycleEdges;

  // Collect the pair of atom indices for every edge of this cycle
  for(const auto& edge : edges) {
    auto edgeFoundPair = boost::edge(
      edge.front(),
      edge.back(),
      graph
    );

    assert(edgeFoundPair.second);

    cycleEdges.insert(edgeFoundPair.first);
  }

  return cycleEdges;
}

unsigned Cycles::numCycleFamilies() const {
  return RDL_getNofURF(_rdlPtr->dataPtr);
}

unsigned Cycles::numCycleFamilies(const AtomIndexType index) const {
  return RDL_getNofURFContainingNode(_rdlPtr->dataPtr, index);
}

//! Returns the number of relevant cycles (RCs)
unsigned Cycles::numRelevantCycles() const {
  return RDL_getNofRC(_rdlPtr->dataPtr);
}

unsigned Cycles::numRelevantCycles(const AtomIndexType index) const {
  return RDL_getNofRCFContainingNode(_rdlPtr->dataPtr, index);
}

Cycles::constIterator Cycles::begin() const {
  return {_rdlPtr};
}

Cycles::constIterator Cycles::end() const {
  return {
    _rdlPtr,
    predicates::All {},
    numRelevantCycles()
  };
}

RDL_data* Cycles::dataPtr() const {
  return _rdlPtr->dataPtr;
}

bool Cycles::operator == (const Cycles& other) const {
  return _rdlPtr == other._rdlPtr;
}

bool Cycles::operator != (const Cycles& other) const {
  return !(*this == other);
}

/* Cycles::RDLDataPtrs */
Cycles::RDLDataPtrs::RDLDataPtrs(const GraphType& sourceGraph) {
  // Initialize a new graph
  graphPtr = RDL_initNewGraph(boost::num_vertices(sourceGraph));

  // Copy vertices (without bond type information)
  auto iterPair = boost::edges(sourceGraph);
  for(auto& iter = iterPair.first; iter != iterPair.second; ++iter) {
    auto edgeAddResult = RDL_addUEdge(
      graphPtr,
      boost::source(*iter, sourceGraph),
      boost::target(*iter, sourceGraph)
    );

    assert(edgeAddResult != RDL_INVALID_RESULT || edgeAddResult != RDL_DUPLICATE_EDGE);
  }

  // Calculate, and assure yourself of success
  dataPtr = RDL_calculate(graphPtr);
  assert(dataPtr != nullptr);
}

Cycles::RDLDataPtrs::~RDLDataPtrs() {
  // Awfully enough, calling this frees both dataPtr and graphPtr
  RDL_deleteData(dataPtr);
}

/* predicates */
bool Cycles::predicates::All::operator() (const RDL_cycle* const) const {
  return true;
}

Cycles::predicates::SizeLessThan::SizeLessThan(unsigned threshold) : threshold {threshold} {}

bool Cycles::predicates::SizeLessThan::operator() (const RDL_cycle* const cyclePtr) const {
  return (cyclePtr -> weight) < threshold;
}

Cycles::predicates::ContainsIndex::ContainsIndex(AtomIndexType soughtIndex) : soughtIndex {soughtIndex} {}

bool Cycles::predicates::ContainsIndex::operator() (const RDL_cycle* const cyclePtr) const {
  for(unsigned i = 0; i < cyclePtr->weight; ++i) {
    if(
        cyclePtr -> edges[i][0] == soughtIndex
        || cyclePtr -> edges[i][1] == soughtIndex
      ) {
      return true;
      break;
    }
  }

  return false;
}

bool Cycles::predicates::ConsistsOf::operator() (const RDL_cycle* const cyclePtr) const {
  if(cyclePtr->weight != indices.size()) {
    return false;
  }

  for(unsigned i = 0; i < cyclePtr->weight; ++i) {
    if(
      !indices.count(cyclePtr->edges[i][0])
      || !indices.count(cyclePtr->edges[i][1])
    ) {
      return false;
    }
  }

  return true;
}

/* Cycles::constIterator */
Cycles::constIterator::constIterator(
  const std::shared_ptr<RDLDataPtrs>& dataPtr,
  Cycles::constIterator::PredicateType cyclePredicate,
  unsigned rCycleIndex
) : _rdlPtr {dataPtr} {
  _cyclePtr = std::make_shared<RDLCyclePtrs>(*dataPtr);
  _cyclePermissiblePredicate = cyclePredicate;

  if(rCycleIndex == 0) { // begin constructor
    // In case first cycle is not permissible, advance to first permissible
    if(
      !RDL_cycleIteratorAtEnd(_cyclePtr->cycleIterPtr)
      && !_cyclePermissiblePredicate(_cyclePtr->cyclePtr)
    ) {
      ++(*this);
    }
  } else { // end constructor, advance to end
    while(_rCycleIndex < rCycleIndex) {
      ++(*this);
    }

    assert(_rCycleIndex == rCycleIndex);
  }
}

Cycles::constIterator& Cycles::constIterator::operator ++ () {
  if(!RDL_cycleIteratorAtEnd(_cyclePtr->cycleIterPtr)) {
    do {
      _cyclePtr->advance();
      ++_rCycleIndex;
    } while(
      !RDL_cycleIteratorAtEnd(_cyclePtr->cycleIterPtr)
      && !_cyclePermissiblePredicate(_cyclePtr->cyclePtr)
    );
  } else {
    throw std::logic_error("Advancing Cycles::constIterator past end");
  }

  return *this;
}

Cycles::constIterator Cycles::constIterator::operator ++ (int) {
  constIterator retval = *this;
  ++(*this);
  return retval;
}

RDL_cycle* Cycles::constIterator::operator * () const {
  if(_cyclePtr->cyclePtr == nullptr) {
    throw std::range_error("Dereferencing Cycles::constIterator end iterator");
  }

  return _cyclePtr->cyclePtr;
}

//! Must be constructed from same Cycles base and at same RC
bool Cycles::constIterator::operator == (const Cycles::constIterator& other) const {
  return (
    _rdlPtr == other._rdlPtr
    && _rCycleIndex == other._rCycleIndex
  );
}

bool Cycles::constIterator::operator != (const Cycles::constIterator& other) const {
  return !(*this == other);
}

/* Cycles::constIterator::RDLCyclePtrs implementation */
Cycles::constIterator::RDLCyclePtrs::RDLCyclePtrs(const RDLDataPtrs& dataPtrs) {
  cycleIterPtr = RDL_getRCyclesIterator(dataPtrs.dataPtr);
  if(RDL_cycleIteratorAtEnd(cycleIterPtr)) {
    cyclePtr = nullptr;
  } else {
    cyclePtr = RDL_cycleIteratorGetCycle(cycleIterPtr);
  }
}

Cycles::constIterator::RDLCyclePtrs::~RDLCyclePtrs() {
  if(cyclePtr != nullptr) {
    RDL_deleteCycle(cyclePtr);
  }

  RDL_deleteCycleIterator(cycleIterPtr);
}

void Cycles::constIterator::RDLCyclePtrs::advance() {
  assert(!RDL_cycleIteratorAtEnd(cycleIterPtr));

  RDL_deleteCycle(cyclePtr);
  cyclePtr = nullptr;

  RDL_cycleIteratorNext(cycleIterPtr);

  if(!RDL_cycleIteratorAtEnd(cycleIterPtr)) {
    cyclePtr = RDL_cycleIteratorGetCycle(cycleIterPtr);
  }
}

std::map<AtomIndexType, unsigned> makeSmallestCycleMap(const Cycles& cycleData) {
  std::map<AtomIndexType, unsigned> smallestCycle;

  for(const auto cyclePtr : cycleData) {
    const unsigned cycleSize = Cycles::size(cyclePtr);

    for(const auto& edgeIndices : Cycles::edgeVertices(cyclePtr)) {
      for(const auto& index: edgeIndices) {
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
  const Cycles::EdgeList& edgeSet
) {
  // copy the edges so we can modify
  auto edgeDescriptors = edgeSet;

  auto firstEdgeIter = edgeDescriptors.begin();
  // Initialize with first edge descriptor's indices
  std::vector<AtomIndexType> indexSequence {
    firstEdgeIter->front(),
    firstEdgeIter->back()
  };

  edgeDescriptors.erase(firstEdgeIter);
  // firstEdgeIter is now invalid!

  while(!edgeDescriptors.empty()) {
    for(
      auto edgeIter = edgeDescriptors.begin();
      edgeIter != edgeDescriptors.end();
      ++edgeIter
    ) {
      if(edgeIter->front() == indexSequence.back()) {
        indexSequence.emplace_back(
          edgeIter->back()
        );
        edgeDescriptors.erase(edgeIter);
        break;
      }

      if(edgeIter->back() == indexSequence.back()) {
        indexSequence.emplace_back(
          edgeIter->front()
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
  const temple::TinySet<GraphType::edge_descriptor>& edgeSet,
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

} // namespace molassembler
