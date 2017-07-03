#ifndef INCLUDE_URF_DATA_H
#define INCLUDE_URF_DATA_H

#include "common_typedefs.h"
#include "RangeForTemporary.h"

// RDL
#include "RingDecomposerLib/RingDecomposerLib.h"

/* ISSUE
 * - Memory management issues -> use-after-free, destructor of one iterator is 
 *   called twice...
 * - Now that RangeForTemporary uses move, the trivial move constructor does a
 *   member-wise copy, on elimination of the temporary the cycleiterator both
 *   are pointing to is deleted
 */

namespace MoleculeManip {

class CycleIterator;

class CycleData {
public:
  friend class CycleIterator;

private:
  const GraphType& _baseGraphReference;
  RDL_graph* _graphPtr;
  RDL_data* _dataPtr;

public:
  // Rule of five methods
  CycleData(const GraphType& sourceGraph);
  CycleData(CycleData&& other);
  CycleData(const CycleData& other) = delete;
  CycleData& operator = (const CycleData& other) = delete;
  ~CycleData();

/* Information */
  //! Returns an iterator proxy object
  CycleIterator getCyclesIterator() const;

  CycleIterator getCyclesIteratorSizeLE(
    const unsigned& maxCycleSize
  ) const;

  CycleIterator getCyclesIteratorContaining(
    const AtomIndexType& containingIndex
  ) const;

  //! Returns the number of unique ring families (URFs)
  unsigned numCycleFamilies() const;

  //! Returns the number of relevant cycles (RCs)
  unsigned numRelevantCycles() const;
};

class CycleIterator {
public:
  friend class CycleData;

private:
  const CycleData& _cycleDataRef;
  const boost::optional<unsigned> _maxCycleSizeOption;
  const boost::optional<AtomIndexType> _containingIndexOption;
  RDL_cycleIterator* _cycleIteratorPtr;

  //! Private constructor so that only CycleData can call it
  CycleIterator(
    const CycleData& cycleData,
    const boost::optional<unsigned>& maxCycleSizeOption = boost::none,
    const boost::optional<AtomIndexType>& containingIndexOption = boost::none
  );

  bool _currentCyclePermissible() const;

public:
  ~CycleIterator();

/* Information */
  bool atEnd() const;
  unsigned cycleSize() const;
  std::set<GraphType::edge_descriptor> getCurrentCycle() const;

/* Modification */
  void advance();
};

// Forward-declare AdjacencyList
// TODO check if this kills the program
class AdjacencyList;

/*!
 * Creates a mapping from atom index to the size of the smallest cycle
 * containing that index. The map does not contain entries for indices not
 * enclosed by a cycle.
 */
std::map<AtomIndexType, unsigned> makeSmallestCycleMap(
  const CycleData& cycleData,
  const AdjacencyList& adjacencies
);

/*!
 * From a set of graph edge descriptors, this function creates one of the two
 * possible vertex index sequences describing the cycle
 */
std::vector<AtomIndexType> makeRingIndexSequence(
  const std::set<GraphType::edge_descriptor>& edgeSet,
  const AdjacencyList& adjacencies
);

/*!
 * Counts the number of planarity enforcing bonds in a set of edge descriptors.
 * Double and aromatic bonds are considered planarity enforcing.
 */
unsigned countPlanarityEnforcingBonds(
  const std::set<GraphType::edge_descriptor>& edgeSet,
  const AdjacencyList& adjacencies
);

} // namespace MoleculeManip

#endif
