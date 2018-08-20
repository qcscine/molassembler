#ifndef INCLUDE_URF_DATA_H
#define INCLUDE_URF_DATA_H

// RDL
#include "RingDecomposerLib.h"

#include "temple/TinySet.h"
#include "temple/Traits.h"

#include "molassembler/detail/RangeForTemporary.h"
#include "molassembler/detail/SharedTypes.h"

/*! @file
 *
 * Contains a wrapper class for the C-style RingDecomposerLib functions so that
 * cycle data can be used in idiomatic C++.
 */

namespace molassembler {

/*! Wrapper class to make working with RDL in C++ more pleasant.
 *
 * Calculated data from a graph is movable, copyable and assignable in all the
 * usual ways, and is therefore suited for caching. Equality comparison of
 * this type and its nested types follows same-base, not equal-base logic of
 * comparison due to management of the C pointer types and the associated
 * allocated memory using shared_ptrs.
 */
class Cycles {
private:
  /*! Safe wrapper around RDL's graph and calculated data pointers
   *
   * Limited operability type to avoid any accidental moves or copies. Manages
   * memory allocation and destruction in all situations.
   */
  struct RDLDataPtrs {
    //! Raw pointers to RDL_graph object
    RDL_graph* graphPtr;

    //! Raw pointer to calculated graph cycle data
    RDL_data* dataPtr;

    RDLDataPtrs() = delete;
    RDLDataPtrs(
      const GraphType& sourceGraph,
      const bool ignoreEtaBonds
    );
    RDLDataPtrs(const RDLDataPtrs& other) = delete;
    RDLDataPtrs(RDLDataPtrs&& other) = delete;
    RDLDataPtrs& operator = (const RDLDataPtrs& other) = delete;
    RDLDataPtrs& operator = (RDLDataPtrs&& other) = delete;
    ~RDLDataPtrs();
  };

  std::shared_ptr<RDLDataPtrs> _rdlPtr;

public:
  Cycles() = default;
  Cycles(const GraphType& sourceGraph, const bool ignoreEtaBonds = false);

  using EdgeList = std::vector<
    std::array<AtomIndexType, 2>
  >;

  //! Returns the size of a cycle (the number of constuting atoms or edges)
  static unsigned size(const RDL_cycle* const cyclePtr);
  //! Returns a list of edges constituted by their edge vertices of a cycle
  static EdgeList edgeVertices(const RDL_cycle* const cyclePtr);
  //! Returns a list of edge descriptors from the original graph of a cycle
  static temple::TinySet<GraphType::edge_descriptor> edges(
    const RDL_cycle* const cyclePtr,
    const GraphType& graph
  );
  //! Returns a list of edge descriptors from an EdgeList
  static temple::TinySet<GraphType::edge_descriptor> edges(
    const EdgeList& edges,
    const GraphType& graph
  );

  //! Returns the number of unique ring families (URFs)
  unsigned numCycleFamilies() const;

  //! Returns the number of unique ring families (URFs) an index is involved in
  unsigned numCycleFamilies(const AtomIndexType index) const;

  //! Returns the number of relevant cycles (RCs)
  unsigned numRelevantCycles() const;

  //! Returns the number of relevant cycles (RCs)
  unsigned numRelevantCycles(const AtomIndexType index) const;

  //! namespace struct
  struct predicates {
    //! Permits all cycles
    struct All {
      bool operator() (const RDL_cycle* const) const;
    };

    //! Limits iterator to cycles smaller than the threshold used in constructor
    struct SizeLessThan {
      const unsigned threshold;

      explicit SizeLessThan(unsigned threshold);

      bool operator() (const RDL_cycle* const cyclePtr) const;
    };

    //! Limits iterator to cycles containing a certain atom index
    struct ContainsIndex {
      AtomIndexType soughtIndex;

      explicit ContainsIndex(AtomIndexType soughtIndex);

      bool operator() (const RDL_cycle* const cyclePtr) const;
    };

    struct ConsistsOf {
      temple::TinySet<AtomIndexType> indices;

      template<typename Container>
      explicit ConsistsOf(const Container& container) {
        static_assert(
          std::is_same<
            temple::traits::getValueType<Container>,
            AtomIndexType
          >::value,
          "ConsistsOf predicate ctor Container must contain AtomIndexType"
        );

        for(const auto index : container) {
          indices.insert(index);
        }
      }

      bool operator() (const RDL_cycle* const cyclePtr) const;
    };
  };

  //! non-modifying forward iterator
  class constIterator {
  public:
    // Iterator traits
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = RDL_cycle*;
    using pointer = value_type;
    using reference = value_type;

    using PredicateType = std::function<bool(const RDL_cycle* const)>;

  private:
    /*! Safe wrapper around RDL's cycle iterator and cycle pointers
     *
     * Limited operability type to avoid any accidental moves or copies. Manages
     * memory allocation and destruction in all situations. Pointer correctness
     * and iterator advancement is also provided.
     */
    struct RDLCyclePtrs {
      RDL_cycleIterator* cycleIterPtr;
      RDL_cycle* cyclePtr;

      RDLCyclePtrs() = delete;
      RDLCyclePtrs(const RDLDataPtrs& dataPtrs);
      RDLCyclePtrs(const RDLCyclePtrs& other) = delete;
      RDLCyclePtrs(RDLCyclePtrs&& other) = delete;
      RDLCyclePtrs& operator = (const RDLCyclePtrs& other) = delete;
      RDLCyclePtrs& operator = (RDLCyclePtrs&& other) = delete;
      ~RDLCyclePtrs();

      /*! Advance internal iterator and cycle state
       *
       * Frees the memory for the current cycle and advances the iterator state.
       * If the iterator is now not at the end of all relevant cycles, then
       * the next cycle is allocated and constructed into cyclePtr. Otherwise,
       * cyclePtr is a nullptr.
       */
      void advance();
    };

    //! Hold an owning reference to the base data to avoid dangling pointers
    std::shared_ptr<RDLDataPtrs> _rdlPtr;
    //! Manage cycle data as shared pointer to permit expected iterator functionality
    std::shared_ptr<RDLCyclePtrs> _cyclePtr;
    //! Holds a predicate function that determines which cycles are permissible
    PredicateType _cyclePermissiblePredicate;
    //! Current position in the full list of relevant cycles, for comparability and construction
    unsigned _rCycleIndex = 0;

  public:
    constIterator() = default;
    constIterator(
      const std::shared_ptr<RDLDataPtrs>& dataPtr,
      PredicateType cyclePredicate = predicates::All {},
      unsigned rCycleIndex = 0
    );

    constIterator& operator ++ ();
    constIterator operator ++ (int);

    RDL_cycle* operator * () const;

    //! Must be constructed from same Cycles base and at same RC to compare equal
    bool operator == (const constIterator& other) const;
    bool operator != (const constIterator& other) const;
  };

  constIterator begin() const;
  constIterator end() const;

  /*! Iterate through cycles satisfying a predicate.
   *
   * PredicateType must be a unary predicate taking a const RDL_cycle* const.
   * Sample predicates are predicates::SizeLessThan and
   * predicates::ContainsIndex.
   */
  template<typename UnaryPredicate>
  RangeForTemporary<constIterator> iterate(UnaryPredicate&& predicate) const {
    return {
      constIterator {_rdlPtr, predicate},
      constIterator {_rdlPtr, predicate, numRelevantCycles()}
    };
  }

  //! Provide access to calculated data
  RDL_data* dataPtr() const;

  //! Must be copy of another to compare equal. Constructed from same graph does not suffice
  bool operator == (const Cycles& other) const;
  bool operator != (const Cycles& other) const;
};

/*!
 * Creates a mapping from atom index to the size of the smallest cycle
 * containing that index. The map does not contain entries for indices not
 * enclosed by a cycle.
 */
std::map<AtomIndexType, unsigned> makeSmallestCycleMap(const Cycles& cycleData);

/*!
 * From a set of graph edge descriptors, this function creates one of the two
 * possible vertex index sequences describing the cycle
 */
std::vector<AtomIndexType> makeRingIndexSequence(
  const Cycles::EdgeList& edgeSet
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
  const temple::TinySet<GraphType::edge_descriptor>& edgeSet,
  const GraphType& graph
);

} // namespace molassembler

#endif
