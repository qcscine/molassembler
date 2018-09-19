// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_CYCLES_H
#define INCLUDE_MOLASSEMBLER_CYCLES_H

// RDL
#include "RingDecomposerLib.h"

#include "temple/TinySet.h"
#include "temple/Traits.h"

#include "molassembler/OuterGraph.h"

#include <functional>
#include <map>

/*! @file
 *
 * @brief Class to explore cyclic structure of molecules
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
public:
  /*! Safe wrapper around RDL's graph and calculated data pointers
   *
   * Limited operability type to avoid any accidental moves or copies. Manages
   * memory allocation and destruction in all situations.
   */
  struct RDLDataPtrs;

  //! Namespacing struct
  struct predicates {
    //! Permits all cycles
    struct All {
      bool operator() (const RDL_cycle* const) const;
    };

    //! Limits iterator to cycles smaller than the threshold used in constructor
    struct SizeLessThan {
      const unsigned threshold;

      explicit SizeLessThan(unsigned passThreshold);

      bool operator() (const RDL_cycle* const cyclePtr) const;
    };

    //! Limits iterator to cycles containing a certain atom index
    struct ContainsIndex {
      AtomIndex soughtIndex;

      explicit ContainsIndex(AtomIndex passSoughtIndex);

      bool operator() (const RDL_cycle* const cyclePtr) const;
    };

    struct ConsistsOf {
      temple::TinySet<AtomIndex> indices;

      template<typename Container>
      explicit ConsistsOf(const Container& container) {
        static_assert(
          std::is_same<
            temple::traits::getValueType<Container>,
            AtomIndex
          >::value,
          "ConsistsOf predicate ctor Container must contain AtomIndex"
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

    /*! Safe wrapper around RDL's cycle iterator and cycle pointers
     *
     * Limited operability type to avoid any accidental moves or copies. Manages
     * memory allocation and destruction in all situations. Pointer correctness
     * and iterator advancement is also provided.
     */
    struct RDLCyclePtrs;

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

  private:
    //! Hold an owning reference to the base data to avoid dangling pointers
    std::shared_ptr<RDLDataPtrs> _rdlPtr;
    //! Manage cycle data as shared pointer to permit expected iterator functionality
    std::shared_ptr<RDLCyclePtrs> _cyclePtr;
    //! Holds a predicate function that determines which cycles are permissible
    PredicateType _cyclePermissiblePredicate;
    //! Current position in the full list of relevant cycles, for comparability and construction
    unsigned _rCycleIndex = 0;
  };

  using EdgeList = std::vector<
    std::array<AtomIndex, 2>
  >;

//!@name Special member functions
//!@{
  Cycles() = default;
  Cycles(const OuterGraph& sourceGraph, bool ignoreEtaBonds = false);
//!@}

//!@name Static member functions
//!@{
  //! Returns the size of a cycle (the number of constuting atoms or edges)
  static unsigned size(const RDL_cycle* const cyclePtr);
  //! Returns a list of edges constituted by their edge vertices of a cycle
  static EdgeList edgeVertices(const RDL_cycle* const cyclePtr);
  //! Returns a list of edge descriptors from the original graph of a cycle
  static temple::TinySet<BondIndex> edges(const RDL_cycle* const cyclePtr);
  //! Returns a list of edge descriptors from an EdgeList
  static temple::TinySet<BondIndex> edges(const EdgeList& edges);
//!@}

//!@name Information
//!@{
  //! Returns the number of unique ring families (URFs)
  unsigned numCycleFamilies() const;

  //! Returns the number of unique ring families (URFs) an index is involved in
  unsigned numCycleFamilies(AtomIndex index) const;

  //! Returns the number of relevant cycles (RCs)
  unsigned numRelevantCycles() const;

  //! Returns the number of relevant cycles (RCs)
  unsigned numRelevantCycles(AtomIndex index) const;

  //! Provide access to calculated data
  RDL_data* dataPtr() const;
//!@}

//!@name Iterators
//!@{
  constIterator begin() const;
  constIterator end() const;

  /*! Provides two iterators that permit iteration through cycles satisfying a predicate.
   *
   * PredicateType must be a unary predicate taking a const RDL_cycle* const.
   * Sample predicates are predicates::SizeLessThan and
   * predicates::ContainsIndex.
   */
  template<typename UnaryPredicate>
  std::pair<constIterator, constIterator> iteratorPair(UnaryPredicate&& predicate) const {
    return {
      constIterator {_rdlPtr, predicate},
      constIterator {_rdlPtr, predicate, numRelevantCycles()}
    };
  }
//!@}

//!@name Operators
//!@{
  //! Must be copy of another to compare equal. Constructed from same graph does not suffice
  bool operator == (const Cycles& other) const;
  bool operator != (const Cycles& other) const;
//!@}

private:
  std::shared_ptr<RDLDataPtrs> _rdlPtr;
};

/*!
 * Creates a mapping from atom index to the size of the smallest cycle
 * containing that index. The map does not contain entries for indices not
 * enclosed by a cycle.
 */
std::map<AtomIndex, unsigned> makeSmallestCycleMap(const Cycles& cycleData);

/*!
 * From a set of graph edge descriptors, this function creates one of the two
 * possible vertex index sequences describing the cycle
 */
std::vector<AtomIndex> makeRingIndexSequence(
  const Cycles::EdgeList& edgeSet
);

std::vector<AtomIndex> centralizeRingIndexSequence(
  std::vector<AtomIndex> ringIndexSequence,
  AtomIndex center
);

/*!
 * Counts the number of planarity enforcing bonds in a set of edge descriptors.
 * Double and aromatic bonds are considered planarity enforcing.
 */
unsigned countPlanarityEnforcingBonds(
  const temple::TinySet<BondIndex>& edgeSet,
  const OuterGraph& graph
);

} // namespace molassembler

#endif
