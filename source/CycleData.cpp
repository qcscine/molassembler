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

void initCycleIter() {
}

/* Destructor */
CycleData::~CycleData() {
  RDL_deleteData(_dataPtr); // this also deletes the graph
}

/* Iterator implementation */
//CycleData::iterator::iterator(
//  const RDL_data* const dataPtr,
//  const unsigned& maxCycleSize,
//  const bool& endPosition
//) : _dataPtr(dataPtr),
//    _maxCycleSize(maxCycleSize),
//    _done(endPosition)
//{
//  // If this is the begin iterator, initialize the cycles iterator
//  if(!endPosition) {
//    _it = RDL_getRCyclesIterator(dataPtr);
//
//    /* If maxCycleSize is not the default, skip ahead to an acceptable cycle
//     * or else all the way to the end
//     */
//    if(maxCycleSize != 0) {
//      bool continueSkipping = true;
//      // Skip ahead to next permissible cycle
//      while(
//       continueSkipping 
//       && RDL_cycleIteratorAtEnd(_it) == 0
//      ) {
//        RDL_cycle* cyclePtr = RDL_cycleIteratorGetCycle(_it);
//        unsigned cycleSize = cyclePtr -> weight;
//        RDL_deleteCycle(cyclePtr);
//
//        continueSkipping = cycleSize > maxCycleSize;
//        if(continueSkipping) {
//          std::cout << "Skipping cycle with size " << cycleSize << std::endl;
//          RDL_cycleIteratorNext(_it);
//        }
//      }
//    }
//  }
//}
//
//CycleData::iterator::iterator(iterator&& other) 
//  : _dataPtr(other._dataPtr),
//    _maxCycleSize(other._maxCycleSize) 
//{
//  std::swap(_it, other._it);
//  std::swap(_done, other._done);
//}
//
//CycleData::iterator::iterator(const iterator& other)
//  : _dataPtr(other._dataPtr),
//    _maxCycleSize(other._maxCycleSize),
//    _done(other._done)
//{
//  // Similar to constructor
//
//}
//
//CycleData::iterator::~iterator() {
//  if(_it != nullptr) {
//    RDL_deleteCycleIterator(_it);
//  }
//}
//
//CycleData::iterator& CycleData::iterator::operator ++ () { 
//  if(
//    !_done 
//    && _it != nullptr
//    && RDL_cycleIteratorAtEnd(_it) == 0
//  ) {
//    RDL_cycleIteratorNext(_it);
//    // We aren't done as long as cycleIteratorAtEnd returns zero
//    _done = RDL_cycleIteratorAtEnd(_it) == 0;
//  }
//
//  return *this; 
//}
//
//bool CycleData::iterator::operator == (const iterator& other) const {
//  if(
//    _done == other._done
//    && _done
//  ) {
//    return true;
//  } 
//
//  return _it == other._it;
//}
//
//bool CycleData::iterator::operator != (const iterator& other) const {
//  return !(
//    *this == other
//  );
//}
//
//typename CycleData::BaseIteratorType::reference CycleData::iterator::operator * () const { 
//  IteratorValueType indices;
//
//  if(!_done) {
//    std::cout << "Accessing cycle..." << std::endl;
//    RDL_cycle* cyclePtr = RDL_cycleIteratorGetCycle(_it);
//
//    for(unsigned i = 0; i < cyclePtr->weight; i++) {
//      indices.insert({
//        cyclePtr->edges[i][0],
//        cyclePtr->edges[i][1]
//      });
//    }
//
//    RDL_deleteCycle(cyclePtr);
//  }
//
//  return indices;
//}*/

/* Information */
unsigned CycleData::numCycleFamilies() const {
  return RDL_getNofURF(_dataPtr);
}

//! Returns the number of relevant cycles (RCs)
unsigned CycleData::numRelevantCycles() const {
  return RDL_getNofRC(_dataPtr);
}

CycleIterator::CycleIterator(
  const CycleData& cycleData,
  const unsigned& maxCycleSize
) : _cycleDataRef(cycleData),
    _maxCycleSize(maxCycleSize) 
{
  // Initialize the cycle iterator
  _cycleIteratorPtr = RDL_getRCyclesIterator(cycleData._dataPtr);

  // Skip to first permissible ring per cycle size restrictions
  if(maxCycleSize != 0) {
    bool continueSkipping = true;
    // Skip ahead to next permissible cycle
    while(
     continueSkipping 
     && RDL_cycleIteratorAtEnd(_cycleIteratorPtr) == 0
    ) {
      RDL_cycle* cyclePtr = RDL_cycleIteratorGetCycle(_cycleIteratorPtr);
      unsigned cycleSize = cyclePtr -> weight;
      RDL_deleteCycle(cyclePtr);

      continueSkipping = cycleSize > maxCycleSize;
      if(continueSkipping) {
        RDL_cycleIteratorNext(_cycleIteratorPtr);
      }
    }
  }
}

bool CycleIterator::atEnd() const {
  return RDL_cycleIteratorAtEnd(_cycleIteratorPtr) != 0;
}

void CycleIterator::advance() {
  assert(!atEnd());
  RDL_cycleIteratorNext(_cycleIteratorPtr);
}

std::vector<GraphType::edge_descriptor> CycleIterator::getCurrentCycle() const {
  std::vector<GraphType::edge_descriptor> cycleEdges;

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

      cycleEdges.push_back(edgeIndexFoundPair.first);
    }

    RDL_deleteCycle(cyclePtr);
  }

  return cycleEdges;
}

CycleIterator::~CycleIterator() {
  RDL_deleteCycleIterator(_cycleIteratorPtr);
}

CycleIterator CycleData::getIterator(const unsigned& maxCycleSize) const {
  return CycleIterator(*this, maxCycleSize);
}

} // namespace MoleculeManip
