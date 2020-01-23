/*! @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Class to explore cyclic structure of molecules
 *
 * Contains a wrapper class for the C-style RingDecomposerLib functions so that
 * cycle data can be used in idiomatic C++.
 */

#ifndef INCLUDE_MOLASSEMBLER_CYCLES_H
#define INCLUDE_MOLASSEMBLER_CYCLES_H

#include "boost/functional/hash.hpp"
#include "boost/optional/optional_fwd.hpp"

#include "molassembler/Types.h"

#include <functional>
#include <unordered_map>

struct RDL_data;

namespace Scine {
namespace molassembler {

// Forward-declarations
class Graph;
class PrivateGraph;

/*!
 * @brief Wrapper class to make working with RDL in C++ more pleasant.
 *
 * Calculated data from a graph is movable, copyable and assignable in all the
 * usual ways, and is therefore suited for caching. Equality comparison of
 * this type and its nested types follows same-base, not equal-base logic of
 * comparison due to management of the C pointer types and the associated
 * allocated memory using shared_ptrs.
 *
 * @note Keep in mind that if a node is part of a URF, that does NOT imply that
 * it is also part of each RC of that URF.
 */
class MASM_EXPORT Cycles {
public:
  /*!
   * @brief Safe wrapper around RDL's graph and calculated data pointers
   *
   * Limited operability type to avoid any accidental moves or copies. Manages
   * memory allocation and destruction in all situations.
   */
  struct RdlDataPtrs;


  /*!
   * @brief Safe wrapper around RDL's cycle iterator and cycle pointers
   *
   * Limited operability type to avoid any accidental moves or copies. Manages
   * memory allocation and destruction in all situations. Pointer correctness
   * and iterator advancement is also provided.
   */
  struct RdlCyclePtrs;

  //! Iterator for all relevant cycles of the graph
  class AllCyclesIterator {
  public:
    // Iterator traits
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = const std::vector<BondIndex>&;
    using pointer = const std::vector<BondIndex>*;
    using reference = value_type;

    AllCyclesIterator() = delete;
    AllCyclesIterator(
      const std::shared_ptr<RdlDataPtrs>& dataPtr,
      unsigned rCycleIndex = 0
    );

    /* Rule of five members */
    AllCyclesIterator(const AllCyclesIterator& other);
    AllCyclesIterator(AllCyclesIterator&& other) noexcept;
    AllCyclesIterator& operator = (const AllCyclesIterator& other);
    AllCyclesIterator& operator = (AllCyclesIterator&& other) noexcept;
    ~AllCyclesIterator();

    AllCyclesIterator& operator ++ ();
    AllCyclesIterator operator ++ (int);
    value_type operator * () const;
    pointer operator -> () const;


    //! Must be constructed from same Cycles base and at same RC to compare equal
    bool operator == (const AllCyclesIterator& other) const;
    bool operator != (const AllCyclesIterator& other) const;

  private:
    //! Hold an owning reference to the base data to avoid dangling pointers
    std::shared_ptr<RdlDataPtrs> _rdlPtr;
    //! Manage cycle data as shared pointer to permit expected iterator functionality
    std::unique_ptr<RdlCyclePtrs> _cyclePtr;
  };

  //! Iterator for cycles of specific universal ring families
  class UrfIdsCycleIterator {
  public:
    using difference_type = unsigned;
    using value_type = const std::vector<BondIndex>&;
    using pointer = const std::vector<BondIndex>*;
    using reference = value_type;
    using iterator_category = std::forward_iterator_tag;

    /* Rule of five members */
    UrfIdsCycleIterator(const UrfIdsCycleIterator& other);
    UrfIdsCycleIterator(UrfIdsCycleIterator&& other) noexcept;
    UrfIdsCycleIterator& operator = (const UrfIdsCycleIterator& other);
    UrfIdsCycleIterator& operator = (UrfIdsCycleIterator&& other) noexcept;
    ~UrfIdsCycleIterator();

    UrfIdsCycleIterator() = delete;

    /* Constructors */
    UrfIdsCycleIterator(
      AtomIndex soughtIndex,
      const std::shared_ptr<RdlDataPtrs>& dataPtr
    );

    UrfIdsCycleIterator(
      const BondIndex& soughtBond,
      const std::vector<unsigned> urfs,
      const std::shared_ptr<RdlDataPtrs>& dataPtr
    );

    UrfIdsCycleIterator(
      const std::vector<BondIndex>& soughtBonds,
      const std::vector<unsigned> urfs,
      const std::shared_ptr<RdlDataPtrs>& dataPtr
    );

    UrfIdsCycleIterator& operator ++ ();
    UrfIdsCycleIterator operator ++ (int);
    value_type operator * () const;
    pointer operator -> () const;

    void advanceToEnd();

    bool operator == (const UrfIdsCycleIterator& other) const;
    bool operator != (const UrfIdsCycleIterator& other) const;

  private:
    struct UrfHelper;

    std::shared_ptr<RdlDataPtrs> _rdlPtr;
    std::unique_ptr<UrfHelper> _urfsPtr;
    std::unique_ptr<RdlCyclePtrs> _cyclePtr;

    void _advanceToNextPermissibleCycle();
    void _initializeCyclesFromURFID();
    void _matchCycleState(const UrfIdsCycleIterator& other);
  };

//!@name Special member functions
//!@{
  /*! @brief Constructor from outer graph
   *
   * @complexity{Approximately linear in the number of bonds in cycles}
   */
  Cycles(const Graph& sourceGraph, bool ignoreEtaBonds = true);
  //! @overload
  Cycles(const PrivateGraph& innerGraph, bool ignoreEtaBonds = true);
//!@}

//!@name Information
//!@{
  /*! @brief Returns the number of unique ring families (URFs)
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned numCycleFamilies() const;

  /*! @brief Returns the number of unique ring families (URFs) an index is involved in
   *
   * @complexity{@math{\Theta(U)} where @math{U} is the number of unique ring
   * families in the molecule}
   */
  unsigned numCycleFamilies(AtomIndex index) const;

  /*! @brief Returns the number of relevant cycles (RCs)
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned numRelevantCycles() const;

  /*! @brief Returns the number of relevant cycles (RCs)
   *
   * @complexity{@math{\Theta(U)} where @math{U} is the number of unique ring
   * families in the molecule}
   */
  unsigned numRelevantCycles(AtomIndex index) const;

  //! Provide access to calculated data
  RDL_data* dataPtr() const;
//!@}

//!@name Iterators
//!@{
  AllCyclesIterator begin() const;
  AllCyclesIterator end() const;
//!@}

//!@name Ranges
//!@{
  /*! @brief Range of relevant cycles containing an atom
   *
   * @complexity{@math{O(U)} where @math{U} is the number of unique ring
   * families of the molecule}
   */
  std::pair<UrfIdsCycleIterator, UrfIdsCycleIterator> containing(AtomIndex atom) const;
  /*! @brief Range of relevant cycles containing a bond
   *
   * @complexity{@math{\Theta(1)}}
   */
  std::pair<UrfIdsCycleIterator, UrfIdsCycleIterator> containing(const BondIndex& bond) const;
  /*! @brief Range of relevant cycles containing several bonds
   *
   * @complexity{@math{\Theta(B)} where @math{B} is the number of bonds in the
   * parameters}
   */
  std::pair<UrfIdsCycleIterator, UrfIdsCycleIterator> containing(const std::vector<BondIndex>& bonds) const;
//!@}

//!@name Operators
//!@{
  //! Must be copy of another to compare equal. Constructed from same graph does not suffice
  bool operator == (const Cycles& other) const;
  bool operator != (const Cycles& other) const;
//!@}

private:
  std::shared_ptr<RdlDataPtrs> _rdlPtr;
  // Map from BondIndex to ordered list of its URF IDs
  std::unordered_map<BondIndex, std::vector<unsigned>, boost::hash<BondIndex>> _urfMap;
};

/*! @brief Yields the size of the smallest cycle containing an atom
 *
 * @complexity{@math{O(U + C)} where @math{U} is the number of unique ring
 * families of the molecule and @math{C} is the number of cycles containing
 * @p atom}
 *
 * @warning Do not use this a lot. Consider makeSmallestCycleMap() instead.
 */
MASM_EXPORT boost::optional<unsigned> smallestCycleContaining(AtomIndex atom, const Cycles& cycles);

/*!
 * @brief Creates a mapping from atom index to the size of the smallest cycle
 *   containing that index.
 *
 * @complexity{@math{\Theta(R)} where @math{R} is the number of relevant cycles
 * in the molecule}
 *
 * @note The map does not contain entries for indices not enclosed by a cycle.
 */
MASM_EXPORT std::unordered_map<AtomIndex, unsigned> makeSmallestCycleMap(const Cycles& cycleData);

/*! @brief Create cycle vertex sequence from unordered edges
 *
 * From a set of unordered graph edge descriptors, this function creates one of
 * the two possible vertex index sequences describing the cycle.
 *
 * @complexity{@math{O(E^2)} worst case}
 */
MASM_EXPORT std::vector<AtomIndex> makeRingIndexSequence(
  std::vector<BondIndex> edgeDescriptors
);

/*! @brief Centralize a cycle vertex sequence at a particular vertex
 *
 * @complexity{@math{O(N)}}
 */
MASM_EXPORT std::vector<AtomIndex> centralizeRingIndexSequence(
  std::vector<AtomIndex> ringIndexSequence,
  AtomIndex center
);

/*! @brief Count the number of planarity enforcing bonds
 *
 * Counts the number of planarity enforcing bonds in a set of edge descriptors.
 * Double bonds are considered planarity enforcing.
 *
 * @complexity{@math{O(N)}}
 */
MASM_EXPORT unsigned countPlanarityEnforcingBonds(
  const std::vector<BondIndex>& edgeSet,
  const Graph& graph
);

} // namespace molassembler
} // namespace Scine

#endif
