#ifndef INCLUDE_URF_DATA_H
#define INCLUDE_URF_DATA_H

#include "common_typedefs.h"
#include "RangeForTemporary.h"
#include "temple/TinySet.h"

// RDL
#include "RingDecomposerLib/RingDecomposerLib.h"

/*! @file
 *
 * Contains a wrapper class for the C-style RingDecomposerLib functions so that
 * cycle data can be used in idiomatic C++.
 */

/* ISSUE
 * - TODO is this still relevant / unfixed?
 * - Memory management issues -> use-after-free, destructor of one iterator is 
 *   called twice...
 * - Now that RangeForTemporary uses move, the trivial move constructor does a
 *   member-wise copy, on elimination of the temporary the cycleiterator both
 *   are pointing to is deleted
 */

namespace molassembler {

// Pre-declare CycleIterator so that it can be friended
class CycleIterator;

/*!
 * RingDecomposerLIb wrapper class so that working with cycle data is not
 * polluted by a slew of un-copyable and odd ANSI C types and functions.
 */
class CycleData {
public:
  friend class CycleIterator;

private:
  //! Refrential access to base molecular graph
  const GraphType& _baseGraphReference;

  //! Raw pointers to RDL_graph object
  RDL_graph* _graphPtr;

  //! Raw pointer to calculated graph cycle data
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

  /*!
   * Returns an iterator proxy object for cycles that are smaller or equal to
   * the passed cycle size
   */
  CycleIterator getCyclesIteratorSizeLE(
    const unsigned& maxCycleSize
  ) const;

  //!  Returns an iterator proxy object for cycles that contain a specific index.
  CycleIterator getCyclesIteratorContaining(
    const AtomIndexType& containingIndex
  ) const;

  //! Returns the number of unique ring families (URFs)
  unsigned numCycleFamilies() const;

  //! Returns the number of unique ring families (URFs) an index is involved in
  unsigned numCycleFamilies(const AtomIndexType& index) const;

  //! Returns the number of relevant cycles (RCs)
  unsigned numRelevantCycles() const;

  //! Returns the number of relevant cycles (RCs)
  unsigned numRelevantCycles(const AtomIndexType& index) const;

  RDL_data* getDataPtr();
};

/*!
 * Iterator-like class for going through cycles in the generated data from the
 * original graph.
 */
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

  using EdgeSet = temple::TinySet<GraphType::edge_descriptor>;

/* Information */
  bool atEnd() const;
  unsigned cycleSize() const;
  EdgeSet getCurrentCycle() const;

/* Modification */
  void advance();
};

// Forward-declare Molecule
class Molecule;

/*!
 * Creates a mapping from atom index to the size of the smallest cycle
 * containing that index. The map does not contain entries for indices not
 * enclosed by a cycle.
 */
std::map<AtomIndexType, unsigned> makeSmallestCycleMap(
  const CycleData& cycleData,
  const GraphType& graph
);

/*!
 * From a set of graph edge descriptors, this function creates one of the two
 * possible vertex index sequences describing the cycle
 */
std::vector<AtomIndexType> makeRingIndexSequence(
  const CycleIterator::EdgeSet& edgeSet,
  const GraphType& graph
);

std::vector<AtomIndexType> centralizeRingIndexSequence(
  std::vector<AtomIndexType> ringIndexSequence,
  const AtomIndexType& center
);

/*!
 * Counts the number of planarity enforcing bonds in a set of edge descriptors.
 * Double and aromatic bonds are considered planarity enforcing.
 */
unsigned countPlanarityEnforcingBonds(
  const CycleIterator::EdgeSet& edgeSet,
  const GraphType& graph
);

} // namespace molassembler

#endif
