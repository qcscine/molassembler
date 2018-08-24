#include "molassembler/Layers/Outer.h"

#include "molassembler/Layers/Bridge.h"

namespace molassembler {

/* Implementation for iterator facade */
template<typename T, bool isVertexInitialized>
OuterGraph::InnerBasedIterator<T, isVertexInitialized>::InnerBasedIterator(
  InnerBasedIterator&& other
) noexcept = default;

template<typename T, bool isVertexInitialized>
OuterGraph::InnerBasedIterator<T, isVertexInitialized>&
OuterGraph::InnerBasedIterator<T, isVertexInitialized>::operator = (InnerBasedIterator&& other) noexcept = default;

template<typename T, bool isVertexInitialized>
OuterGraph::InnerBasedIterator<T, isVertexInitialized>::InnerBasedIterator(
  const InnerBasedIterator& other
) : _pImpl (
  std::make_unique<Impl>(*other._pImpl)
) { }

template<typename T, bool isVertexInitialized>
OuterGraph::InnerBasedIterator<T, isVertexInitialized>&
OuterGraph::InnerBasedIterator<T, isVertexInitialized>::operator = (const InnerBasedIterator& other) {
  *_pImpl = *other._pImpl;
}

template<typename T, bool isVertexInitialized>
OuterGraph::InnerBasedIterator<T, isVertexInitialized>::~InnerBasedIterator() = default;

template<typename T, bool isVertexInitialized>
template<bool Dependent, std::enable_if_t<!Dependent, int>...>
OuterGraph::InnerBasedIterator<T, isVertexInitialized>::InnerBasedIterator(
  const InnerGraph& inner,
  bool begin
) : _pImpl (
  std::make_unique<Impl>(inner, begin)
) {}

template<typename T, bool isVertexInitialized>
template<bool Dependent, std::enable_if_t<Dependent, int>...>
OuterGraph::InnerBasedIterator<T, isVertexInitialized>::InnerBasedIterator(
  const AtomIndex a,
  const InnerGraph& inner,
  bool begin
) : _pImpl (
  std::make_unique<Impl>(a, inner, begin)
) {}

template<typename T, bool isVertexInitialized>
OuterGraph::InnerBasedIterator<T, isVertexInitialized>& OuterGraph::InnerBasedIterator<T, isVertexInitialized>::operator ++ () {
  return ++(*_pImpl);
}

template<typename T, bool isVertexInitialized>
OuterGraph::InnerBasedIterator<T, isVertexInitialized> OuterGraph::InnerBasedIterator<T, isVertexInitialized>::operator ++ (int) {
  auto copy = *this;
  ++(*_pImpl);
  return copy;
}

template<typename T, bool isVertexInitialized>
typename OuterGraph::InnerBasedIterator<T, isVertexInitialized>::value_type OuterGraph::InnerBasedIterator<T, isVertexInitialized>::operator * () const {
  return *(*_pImpl);
}

template<typename T, bool isVertexInitialized>
typename OuterGraph::InnerBasedIterator<T, isVertexInitialized>::difference_type OuterGraph::InnerBasedIterator<T, isVertexInitialized>::operator - (const InnerBasedIterator<T, isVertexInitialized>& other) const {
  return *_pImpl - *other._pImpl;
}

template<typename T, bool isVertexInitialized>
bool OuterGraph::InnerBasedIterator<T, isVertexInitialized>::operator == (const InnerBasedIterator<T, isVertexInitialized>& other) const {
  return *_pImpl == *other._pImpl;
}

template<typename T, bool isVertexInitialized>
bool OuterGraph::InnerBasedIterator<T, isVertexInitialized>::operator != (const InnerBasedIterator<T, isVertexInitialized>& other) const {
  return !(*_pImpl == *other._pImpl);
}

/* Impl base class for quicker implementation of each iterator Impl */
template<typename Iterator>
struct BaseIteratorWrapper {
  Iterator iterator;

  BaseIteratorWrapper() = default;
  BaseIteratorWrapper(Iterator passIterator) : iterator(std::move(passIterator)) {}

  BaseIteratorWrapper& operator ++ () {
    ++iterator;
    return *this;
  }

  /* BaseIteratorWrapper operator ++ (int) is unneeded since a pImpled iterator
   * wrapper never calls it
   */

  typename Iterator::difference_type operator - (const BaseIteratorWrapper& other) const {
    return iterator - other.iterator;
  }

  bool operator == (const BaseIteratorWrapper& other) const {
    return iterator == other.iterator;
  }
};


/* Impl definitions and missing functions */
template<>
struct OuterGraph::AtomIterator::Impl
  : public BaseIteratorWrapper<InnerGraph::BGLType::vertex_iterator>
{
  Impl(const InnerGraph& inner, bool begin) {
    auto vertexIterators = inner.vertices();
    if(begin) {
      iterator = std::move(vertexIterators.first);
    } else {
      iterator = std::move(vertexIterators.second);
    }
  }

  AtomIndex operator * () const {
    // OuterGraph::AtomIndex and InnerGraph::Vertex are the same type
    return *iterator;
  }
};

template<>
struct OuterGraph::BondIterator::Impl
  : public BaseIteratorWrapper<InnerGraph::BGLType::edge_iterator>
{
  const InnerGraph* innerPtr;

  Impl(const InnerGraph& inner, bool begin) : innerPtr(&inner) {
    auto edgeIterators = inner.edges();
    if(begin) {
      iterator = std::move(edgeIterators.first);
    } else {
      iterator = std::move(edgeIterators.second);
    }
  }

  BondIndex operator * () const {
    return toOuter(*iterator, *innerPtr);
  }
};

template<>
struct OuterGraph::AdjacencyIterator::Impl
  : public BaseIteratorWrapper<InnerGraph::BGLType::adjacency_iterator>
{
  Impl(
    const AtomIndex a,
    const InnerGraph& inner,
    bool begin
  ) {
    auto adjacencyIterators = inner.adjacents(a);
    if(begin) {
      iterator = std::move(adjacencyIterators.first);
    } else {
      iterator = std::move(adjacencyIterators.second);
    }
  }

  AtomIndex operator * () const {
    return *iterator;
  }
};

template<>
struct OuterGraph::IncidentEdgesIterator::Impl
  : public BaseIteratorWrapper<InnerGraph::BGLType::out_edge_iterator>
{
  const InnerGraph* innerPtr;

  Impl(const AtomIndex a, const InnerGraph& inner, bool begin) : innerPtr(&inner) {
    auto edgeIterators = inner.edges(a);
    if(begin) {
      iterator = std::move(edgeIterators.first);
    } else {
      iterator = std::move(edgeIterators.second);
    }
  }

  BondIndex operator * () const {
    return toOuter(*iterator, *innerPtr);
  }
};

} // namespace molassembler
