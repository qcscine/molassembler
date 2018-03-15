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

namespace molassembler {

class CycleData {
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

  // Iterator declaration
  using IteratorValueType = std::vector<GraphType::edge_descriptor>;

  using BaseIteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    IteratorValueType,                 // value_type
    unsigned,                  // difference_type
    const IteratorValueType*,          // pointer
    IteratorValueType                  // reference
  >;

  class iterator : public BaseIteratorType {
  private:
    const CycleData& _cycleDataRef;
    // Constant state
    const boost::optional<unsigned> _maxCycleSizeOption;
    const boost::optional<AtomIndexType>& _containingIndexOption;

    // State
    RDL_cycleIterator* _cycleIteratorPtr;
    unsigned _index = 0;

    bool _currentCyclePermissible() const;

    inline bool _atEnd() const;

  public:
    iterator(
      const CycleData& cycleData,
      const boost::optional<unsigned>& maxCycleSizeOption,
      const boost::optional<AtomIndexType>& containingIndexOption,
      const bool& endIterator = false
    );

    // Copy constructor
    iterator(const iterator& other);

    // Move constructor is also a deep copy
    iterator(iterator&& other);

    // Copy assignment
    iterator& operator = (const iterator& other) = delete;

    ~iterator();

    iterator& operator ++ ();

    iterator operator++ (int);

    bool operator == (iterator other) const;
    
    bool operator != (iterator other) const;

    typename BaseIteratorType::value_type operator * () const;
  };

/* Information */
  //! Returns an iterator proxy object
  RangeForTemporary<iterator> iterateCyclesSmallerThan(
    const unsigned& smallerSize
  ) const;
  RangeForTemporary<iterator> iterateCyclesContaining(
    const AtomIndexType& containingIndex
  ) const;

  iterator begin() const;
  iterator end() const;

  //! Returns the number of unique ring families (URFs)
  unsigned numCycleFamilies() const;

  //! Returns the number of relevant cycles (RCs)
  unsigned numRelevantCycles() const;
};

} // namespace molassembler

#endif
