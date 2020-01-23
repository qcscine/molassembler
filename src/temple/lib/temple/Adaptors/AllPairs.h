/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides pair-generation within a single container or two.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ALL_PAIRS_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ALL_PAIRS_ADAPTOR_H

#include "temple/ContainerTraits.h"
#include "temple/Binding.h"

#include <tuple>

namespace Scine {
namespace temple {
namespace adaptors {
namespace detail {

template<class Base>
struct EmptySizeSupplier {};

template<class Base>
struct SingleContainerPairsSizeSupplier {
  std::size_t size() const {
    const std::size_t N = static_cast<const Base&>(*this).container.size();
    return N * (N - 1) / 2;
  }
};

template<class Base>
struct TwoContainersPairsSizeSupplier {
  std::size_t size() const {
    const auto& base = static_cast<const Base&>(*this);
    return base.containerT.size() * base.containerU.size();
  }
};

template<class Container>
struct SingleContainerPairsGenerator
  : public std::conditional_t<
    traits::hasSize<Container>::value,
    SingleContainerPairsSizeSupplier<SingleContainerPairsGenerator<Container>>,
    EmptySizeSupplier<SingleContainerPairsGenerator<Container>>
  >
{
//!@name Types
//!@{
  // See tricks documentation
  using BoundContainer = typename Binding<Container>::type;

  using ContainerValueType = decltype(
    *std::begin(
      std::declval<const Container>()
    )
  );

  using PairType = std::pair<ContainerValueType, ContainerValueType>;

  using ContainerIteratorType = decltype(
    std::begin(std::declval<const Container>())
  );
//!@}

//!@name Public members
//!@{
  BoundContainer container;
//!@}

//!@name Special member functions
//!@{
  explicit SingleContainerPairsGenerator(Container&& passContainer)
    : container(std::forward<Container>(passContainer)) {}
//!@}

//!@name Iterators
//!@{
  template<class ContainerIterator>
  class iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = PairType;
    using difference_type = int;
    using pointer = const PairType*;
    using reference = const PairType&;

    iterator() = default;
    iterator(ContainerIterator left, ContainerIterator right, ContainerIterator end)
      : _left(std::move(left)),
        _right(std::move(right)),
        _end(std::move(end))
    {}

    // Prefix increment
    iterator& operator ++ () {
      ++_right;
      if(_right == _end) {
        ++_left;
        _right = _left;
        ++_right;
      }

      return *this;
    }

    // Postfix increment
    iterator operator ++ (int) {
      iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const iterator& other) const {
      return _left == other._left && _right == other._right;
    }

    bool operator != (const iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {*_left, *_right};
    }

  private:
    ContainerIterator _left, _right, _end;
  };

  iterator<ContainerIteratorType> begin() const {
    auto maybeNextToBegin = std::begin(container);
    if(maybeNextToBegin != std::end(container)) {
      ++maybeNextToBegin;
    }

    return {
      std::begin(container),
      std::move(maybeNextToBegin),
      std::end(container)
    };
  }

  iterator<ContainerIteratorType> end() const {
    auto maybePriorToEnd = std::end(container);
    if(maybePriorToEnd != std::begin(container)) {
      --maybePriorToEnd;
    }

    return {
      std::move(maybePriorToEnd),
      std::end(container),
      std::end(container)
    };
  }
//!@}
};

template<class ContainerT, class ContainerU>
struct TwoContainersAllPairsGenerator
  : public std::conditional_t<
    traits::hasSize<ContainerT>::value && traits::hasSize<ContainerU>::value,
    TwoContainersPairsSizeSupplier<TwoContainersAllPairsGenerator<ContainerT, ContainerU>>,
    SingleContainerPairsSizeSupplier<TwoContainersAllPairsGenerator<ContainerT, ContainerU>>
  >
{
//!@name Types
//!@{
  using BoundContainerT = typename Binding<ContainerT>::type;
  using BoundContainerU = typename Binding<ContainerU>::type;

  using T = decltype(
    *std::begin(
      std::declval<ContainerT>()
    )
  );

  using U = decltype(
    *std::begin(
      std::declval<ContainerU>()
    )
  );

  using PairType = std::pair<T, U>;

  using ContainerTIterator = decltype(
    std::begin(std::declval<const ContainerT>())
  );
  using ContainerUIterator = decltype(
    std::begin(std::declval<const ContainerU>())
  );
//!@}

//!@name Public members
//!@{
  BoundContainerT containerT;
  BoundContainerU containerU;
//!@}

//!@name Special member functions
//!@{
  TwoContainersAllPairsGenerator(
    ContainerT&& t,
    ContainerU&& u
  ) : containerT(std::forward<ContainerT>(t)),
      containerU(std::forward<ContainerU>(u))
  {}
//!@}

//!@name Iterators
//!@{
  template<class TIterator, class UIterator>
  class iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = PairType;
    using difference_type = int;
    using pointer = const PairType*;
    using reference = const PairType&;

    iterator() = default;
    iterator(
      TIterator tBegin,
      TIterator tEnd,
      UIterator uBegin,
      UIterator uEnd
    ) : _tIter(std::move(tBegin)),
        _tEnd(std::move(tEnd)),
        _uBegin(std::move(uBegin)),
        _uIter(_uBegin),
        _uEnd(std::move(uEnd))
    {}

    iterator& operator ++ () {
      ++_uIter;
      if(_uIter == _uEnd) {
        ++_tIter;
        _uIter = _uBegin;
      }

      return *this;
    }

    iterator operator ++ (int) {
      iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const iterator& other) const {
      return (
        std::tie(_tIter, _tEnd, _uIter, _uEnd)
        == std::tie(other._tIter, other._tEnd, other._uIter, other._uEnd)
      );
    }

    bool operator != (const iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {*_tIter, *_uIter};
    }

  private:
    TIterator _tIter, _tEnd;
    UIterator _uBegin, _uIter, _uEnd;
  };

  iterator<ContainerTIterator, ContainerUIterator> begin() const {
    return {
      std::begin(containerT),
      std::end(containerT),
      std::begin(containerU),
      std::end(containerU)
    };
  }

  iterator<ContainerTIterator, ContainerUIterator> end() const {
    return {
      std::end(containerT),
      std::end(containerT),
      std::begin(containerU),
      std::end(containerU)
    };
  }
//!@}
};

} // namespace detail

template<class Container>
auto allPairs(Container&& container) {
  return detail::SingleContainerPairsGenerator<Container>(
    std::forward<Container>(container)
  );
}

template<class ContainerT, class ContainerU>
auto allPairs(ContainerT&& t, ContainerU&& u) {
  return detail::TwoContainersAllPairsGenerator<ContainerT, ContainerU>(
    std::forward<ContainerT>(t),
    std::forward<ContainerU>(u)
  );
}

} // namespace adaptors
} // namespace temple
} // namespace Scine

#endif
