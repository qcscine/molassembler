#include "CycleData.h"

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

//! Returns the number of relevant cycles (RCs)
unsigned CycleData::numRelevantCycles() const {
  return RDL_getNofRC(_dataPtr);
}


bool CycleData::iterator::_currentCyclePermissible() const {
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

inline bool CycleData::iterator::_atEnd() const {
  return RDL_cycleIteratorAtEnd(_cycleIteratorPtr) != 0;
}

CycleData::iterator::iterator(
  const CycleData& cycleData,
  const boost::optional<unsigned>& maxCycleSizeOption,
  const boost::optional<AtomIndexType>& containingIndexOption,
  const bool& endIterator
) : _cycleDataRef(cycleData),
    _maxCycleSizeOption(maxCycleSizeOption),
    _containingIndexOption(containingIndexOption)
{
  // Initialize the cycle iterator
  _cycleIteratorPtr = RDL_getRCyclesIterator(cycleData._dataPtr);

  if(!endIterator) {
    // Skip to first permissible ring per cycle size restrictions
    while(
      !_atEnd()
      && !_currentCyclePermissible() 
    ) {
      RDL_cycleIteratorNext(_cycleIteratorPtr);
      _index += 1;
    }
  } else {
    // Increment to end
    while(!_atEnd()) {
      RDL_cycleIteratorNext(_cycleIteratorPtr);
      _index += 1;
    }
  }
}

// Copy constructor
CycleData::iterator::iterator(const CycleData::iterator& other)
  : iterator( // delegate to default ctor
    other._cycleDataRef,
    other._maxCycleSizeOption,
    other._containingIndexOption
  )
{
  while(_index < other._index && !_atEnd()) {
    RDL_cycleIteratorNext(_cycleIteratorPtr);
    _index += 1;
  }

  assert(_index == other._index);
}

// Move constructor is also a deep copy
CycleData::iterator::iterator(CycleData::iterator&& other)
  : iterator( // delegate to default ctor
    other._cycleDataRef,
    other._maxCycleSizeOption,
    other._containingIndexOption
  )
{
  while(_index < other._index && !_atEnd()) {
    RDL_cycleIteratorNext(_cycleIteratorPtr);
    _index += 1;
  }

  assert(_index == other._index);
}

CycleData::iterator::~iterator() {
  RDL_deleteCycleIterator(_cycleIteratorPtr);
}

CycleData::iterator& CycleData::iterator::operator ++ () {
  assert(!_atEnd());
  // Advance to next permissible cycle or end
  do {
    RDL_cycleIteratorNext(_cycleIteratorPtr);
    _index += 1;
  } while(
    !_atEnd()
    && !_currentCyclePermissible() 
  );
  return *this; 
}

CycleData::iterator CycleData::iterator::operator++ (int) { 
  iterator retval {*this};
  ++(*this);
  return retval; 
}

bool CycleData::iterator::operator == (iterator other) const {
  return (
    &_cycleDataRef == &other._cycleDataRef // address of cycle is same
    && _maxCycleSizeOption == other._maxCycleSizeOption
    && _containingIndexOption == other._containingIndexOption
    && _index == other._index
  );
}

bool CycleData::iterator::operator != (iterator other) const {
  return (
    &_cycleDataRef != &other._cycleDataRef // address of cycle is unequal
    || _maxCycleSizeOption != other._maxCycleSizeOption
    || _containingIndexOption != other._containingIndexOption
    || _index == other._index
  );
}

typename CycleData::BaseIteratorType::value_type CycleData::iterator::operator * () const { 
  std::vector<GraphType::edge_descriptor> cycleEdges;

  if(!_atEnd()) {
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

      cycleEdges.push_back(edgeIndexFoundPair.first);
    }

    RDL_deleteCycle(cyclePtr);
  }

  return cycleEdges;
}

RangeForTemporary<CycleData::iterator> CycleData::iterateCyclesSmallerThan(
  const unsigned& smallerSize
) const {
  return RangeForTemporary<iterator>(
    iterator(
      *this,
      smallerSize - 1,
      boost::none
    ),
    iterator(
      *this,
      smallerSize - 1,
      boost::none,
      true
    )
  );
}

RangeForTemporary<CycleData::iterator> CycleData::iterateCyclesContaining(
  const AtomIndexType& containingIndex
) const {
  return RangeForTemporary<iterator>(
    iterator(
      *this,
      boost::none,
      containingIndex
    ),
    iterator(
      *this,
      boost::none,
      containingIndex,
      true
    )
  );
}

CycleData::iterator CycleData::begin() const {
  return {
    *this,
    boost::none,
    boost::none
  };
}

CycleData::iterator CycleData::end() const {
  return {
    *this,
    boost::none,
    boost::none,
    true
  };
}

} // namespace MoleculeManip
