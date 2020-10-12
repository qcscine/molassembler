/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides pair-generation within a single container or two.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ALL_PAIRS_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ALL_PAIRS_ADAPTOR_H

#include "Molassembler/Temple/ContainerTraits.h"
#include "Molassembler/Temple/Binding.h"

#include <tuple>

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Adaptors {
namespace Detail {

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
    Traits::hasSize<Container>::value,
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
      : left_(std::move(left)),
        right_(std::move(right)),
        end_(std::move(end))
    {}

    // Prefix increment
    iterator& operator ++ () {
      ++right_;
      if(right_ == end_) {
        ++left_;
        right_ = left_;
        ++right_;
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
      return left_ == other.left_ && right_ == other.right_;
    }

    bool operator != (const iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {*left_, *right_};
    }

  private:
    ContainerIterator left_, right_, end_;
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
    Traits::hasSize<ContainerT>::value && Traits::hasSize<ContainerU>::value,
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
    ) : tIter_(std::move(tBegin)),
        tEnd_(std::move(tEnd)),
        uBegin_(std::move(uBegin)),
        uIter_(uBegin_),
        uEnd_(std::move(uEnd))
    {}

    iterator& operator ++ () {
      ++uIter_;
      if(uIter_ == uEnd_) {
        ++tIter_;
        uIter_ = uBegin_;
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
        std::tie(tIter_, tEnd_, uIter_, uEnd_)
        == std::tie(other.tIter_, other.tEnd_, other.uIter_, other.uEnd_)
      );
    }

    bool operator != (const iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {*tIter_, *uIter_};
    }

  private:
    TIterator tIter_, tEnd_;
    UIterator uBegin_, uIter_, uEnd_;
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

} // namespace Detail

template<class Container>
auto allPairs(Container&& container) {
  return Detail::SingleContainerPairsGenerator<Container>(
    std::forward<Container>(container)
  );
}

template<class ContainerT, class ContainerU>
auto allPairs(ContainerT&& t, ContainerU&& u) {
  return Detail::TwoContainersAllPairsGenerator<ContainerT, ContainerU>(
    std::forward<ContainerT>(t),
    std::forward<ContainerU>(u)
  );
}

} // namespace Adaptors
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
