#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_SEQUENTIAL_PAIRS_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_SEQUENTIAL_PAIRS_ADAPTOR_H

#include "temple/ContainerTraits.h"

/*!@file
 *
 * @brief Provides sequential pair generation from a single container
 */

namespace temple {

namespace adaptors {

namespace detail {

template<class Container>
struct SequentialPairGenerator {
//!@name Types
//!@{
  // See tricks documentation
  using BoundContainer = std::conditional_t<
    std::is_rvalue_reference<Container&&>::value,
    std::decay_t<Container>,
    const Container&
  >;

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

//!@name Member variables
//!@{
  BoundContainer container;
//!@}

//!@name Information
//!@{
  std::enable_if_t<
    traits::hasSize<Container>::value,
    std::size_t
  > size() const {
    if(container.size() > 0) {
      return container.size() - 1;
    }

    return 0;
  }
//!@}

//!@name Special member functions
//!@{
  SequentialPairGenerator(Container&& passContainer)
    : container(std::forward<Container>(passContainer)) {}
//!@}

//!@name Iterators
//!@{
  class Iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = PairType;
    using difference_type = int;
    using pointer = const PairType*;
    using reference = const PairType&;

    Iterator() = default;
    Iterator(ContainerIteratorType left, ContainerIteratorType right)
      : _left(std::move(left)),
        _right(std::move(right))
    {}

    Iterator& operator ++ () {
      ++_left;
      ++_right;
      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return _left == other._left && _right == other._right;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {
        *_left,
        *_right
      };
    }

  private:
    ContainerIteratorType _left, _right;
  };

  Iterator begin() const {
    auto maybeNextToBegin = std::begin(container);
    if(maybeNextToBegin != std::end(container)) {
      ++maybeNextToBegin;
    }

    return {
      std::begin(container),
      std::move(maybeNextToBegin)
    };
  }

  Iterator end() const {
    auto maybePriorToEnd = std::end(container);
    if(maybePriorToEnd != std::begin(container)) {
      --maybePriorToEnd;
    }

    return {
      std::move(maybePriorToEnd),
      std::end(container)
    };
  }
//!@}
};

} // namespace detail

template<class Container>
auto sequentialPairs(Container&& container) {
  return detail::SequentialPairGenerator<Container>(
    std::forward<Container>(container)
  );
}

} // namespace adaptors

} // namespace temple

#endif
