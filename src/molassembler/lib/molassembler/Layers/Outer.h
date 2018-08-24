#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_H

#include "boost/optional/optional_fwd.hpp"
#include "Delib/ElementTypes.h"

#include "molassembler/Layers/Types.h"

#include <memory>

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

namespace molassembler {

// Forward-declare InnerGraph
class InnerGraph;

class OuterGraph {
public:
//!@name Member types
//!@{
  using AtomIndex = std::size_t;

  struct BondIndex {
    AtomIndex first, second;

    BondIndex();
    BondIndex(AtomIndex a, AtomIndex b) : first(a), second(b) {
      if(b < a) {
        std::swap(a, b);
      }
    }

    bool operator < (const BondIndex& other) const;
    bool operator == (const BondIndex& other) const;
  };

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

    // Enable construction from Graph and bool if not vertex initialized
    template<bool Dependent = isVertexInitialized, std::enable_if_t<!Dependent, int>...>
    InnerBasedIterator(const InnerGraph& inner, bool begin);

    // Enable compound construction if vertex initialized
    template<bool Dependent = isVertexInitialized, std::enable_if_t<Dependent, int>...>
    InnerBasedIterator(AtomIndex a, const InnerGraph& inner, bool begin);

    InnerBasedIterator& operator ++ ();
    InnerBasedIterator operator ++ (int);
    value_type operator * () const;
    difference_type operator - (const InnerBasedIterator& other) const;

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

//!@name Information
//!@{
  bool adjacent(AtomIndex a, AtomIndex b) const;
  boost::optional<BondIndex> bond(AtomIndex a, AtomIndex b) const;
  BondType bondType(const BondIndex& a) const;
  bool canRemove(AtomIndex a) const;
  bool canRemove(BondIndex a) const;
  Delib::ElementType elementType(AtomIndex a) const;
  AtomIndex N() const;
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
#else
    std::unique_ptr<InnerGraph> _innerPtr;
#endif
};

class Molecule {
public:
  OuterGraph& graph();

private:
  OuterGraph _graph;
};

} // namespace molassembler

#endif
