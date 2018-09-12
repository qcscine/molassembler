#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_H

#include "boost/optional/optional_fwd.hpp"
#include "Delib/ElementTypes.h"

#include "molassembler/Types.h"

#include <memory>

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

/*!@file
 *
 * @brief Interface class for the molecular graph
 */

// Forward-declarations
namespace Delib {
  class ElementTypeCollection;
} // namespace Delib

namespace molassembler {

// Forward-declare InnerGraph
class InnerGraph;
class Cycles;

class OuterGraph {
public:
//!@name Member types
//!@{
  template<typename T, bool isVertexInitialized>
  class InnerBasedIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = T;
    using reference = T;

    InnerBasedIterator(InnerBasedIterator&& other) noexcept;
    InnerBasedIterator& operator = (InnerBasedIterator&& other) noexcept;
    InnerBasedIterator(const InnerBasedIterator& other);
    InnerBasedIterator& operator = (const InnerBasedIterator& other);
    ~InnerBasedIterator();

    // Default constructor
    InnerBasedIterator();

    // Enable construction from Graph and bool if not vertex initialized
    template<bool Dependent = isVertexInitialized, std::enable_if_t<!Dependent, int>...>
    InnerBasedIterator(const InnerGraph& inner, bool begin);

    // Enable compound construction if vertex initialized
    template<bool Dependent = isVertexInitialized, std::enable_if_t<Dependent, int>...>
    InnerBasedIterator(AtomIndex a, const InnerGraph& inner, bool begin);

    InnerBasedIterator& operator ++ ();
    InnerBasedIterator operator ++ (int);
    value_type operator * () const;

    bool operator == (const InnerBasedIterator& other) const;
    bool operator != (const InnerBasedIterator& other) const;

  private:
    struct Impl;

#ifdef MOLASSEMBLER_ENABLE_PROPAGATE_CONST
    std::experimental::propagate_const<
      std::unique_ptr<Impl>
    > _pImpl;
#else
    std::unique_ptr<Impl> _pImpl;
#endif
  };

  using AtomIterator = InnerBasedIterator<AtomIndex, false>;
  using BondIterator = InnerBasedIterator<BondIndex, false>;
  using AdjacencyIterator = InnerBasedIterator<AtomIndex, true>;
  using IncidentEdgesIterator = InnerBasedIterator<BondIndex, true>;

  template<typename Iter>
  using Range = std::pair<Iter, Iter>;
//!@}

//!@name Special member functions
//!@{
  OuterGraph(OuterGraph&& other) noexcept;
  OuterGraph& operator = (OuterGraph&& other) noexcept;
  OuterGraph(const OuterGraph& other);
  OuterGraph& operator = (const OuterGraph& other);
  ~OuterGraph();
//!@{

//!@name Constructors
//!@{
  OuterGraph();
  //! Wrapping constructor
  explicit OuterGraph(InnerGraph&& inner);
//!@}

//!@name Information
//!@{
  bool adjacent(AtomIndex a, AtomIndex b) const;
  boost::optional<BondIndex> bond(AtomIndex a, AtomIndex b) const;
  BondType bondType(const BondIndex& edge) const;
  bool canRemove(AtomIndex a) const;
  bool canRemove(const BondIndex& edge) const;
  const Cycles& cycles() const;
  unsigned degree(AtomIndex a) const;
  Delib::ElementTypeCollection elementCollection() const;
  Delib::ElementType elementType(AtomIndex a) const;

  AtomIndex N() const;
  unsigned B() const;
//!@}

//!@name Ranges
//!@{
  Range<AtomIterator> atoms() const;
  Range<BondIterator> bonds() const;
  Range<AdjacencyIterator> adjacents(AtomIndex a) const;
  Range<IncidentEdgesIterator> bonds(AtomIndex a) const;
//!@}

  // Access to BGL-level descriptors (for inner layer users only)
  InnerGraph& inner() { return *_innerPtr; }
  const InnerGraph& inner() const { return *_innerPtr; }

private:
#ifdef MOLASSEMBLER_ENABLE_PROPAGATE_CONST
  std::experimental::propagate_const<
    std::unique_ptr<InnerGraph>
  > _innerPtr;

  mutable std::experimental::propagate_const<
    std::unique_ptr<Cycles>
  > _cachedCycles;
#else
  std::unique_ptr<InnerGraph> _innerPtr;
  mutable std::unique_ptr<Cycles> _cachedCycles;
#endif
};

} // namespace molassembler

#endif
