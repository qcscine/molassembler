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
  CycleIterator getIterator(const unsigned& maxCycleSize = 0) const;

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
  const unsigned _maxCycleSize;
  RDL_cycleIterator* _cycleIteratorPtr;

  //! Private constructor so that only CycleData can call it
  CycleIterator(
    const CycleData& cycleData,
    const unsigned& maxCycleSize = 0
  );

public:
  ~CycleIterator();

/* Information */
  bool atEnd() const;
  std::vector<GraphType::edge_descriptor> getCurrentCycle() const;

/* Modification */
  void advance();
};

} // namespace MoleculeManip

#endif
