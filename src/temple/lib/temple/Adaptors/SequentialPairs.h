/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides sequential pair generation from a single container
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_SEQUENTIAL_PAIRS_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_SEQUENTIAL_PAIRS_ADAPTOR_H

#include "temple/ContainerTraits.h"
#include "temple/Binding.h"

#include <vector>

namespace Scine {
namespace temple {
namespace adaptors {
namespace detail {

template<class Container>
struct SequentialPairGenerator : public Binding<Container> {
//!@name Types
//!@{
  using ContainerBinding = Binding<Container>;

  using ContainerIteratorType = decltype(
    std::begin(std::declval<std::add_const_t<std::decay_t<Container>>>())
  );

  using ContainerValueType = decltype(
    *std::declval<ContainerIteratorType>()
  );

  using PairType = std::pair<ContainerValueType, ContainerValueType>;
//!@}

//!@name Information
//!@{
  std::enable_if_t<
    traits::hasSize<Container>::value,
    std::size_t
  > size() const {
    if(ContainerBinding::value.size() > 0) {
      return ContainerBinding::value.size() - 1;
    }

    return 0;
  }
//!@}

//!@name Constructor
//!@{
  using ContainerBinding::ContainerBinding;
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
      : left_(std::move(left)),
        right_(std::move(right))
    {}

    Iterator& operator ++ () {
      ++left_;
      ++right_;
      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return left_ == other.left_ && right_ == other.right_;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {*left_, *right_};
    }

  private:
    ContainerIteratorType left_, right_;
  };

  Iterator begin() const {
    const auto& container = ContainerBinding::value;

    auto maybeNextToBegin = std::begin(container);
    if(maybeNextToBegin != std::end(container)) {
      ++maybeNextToBegin;
    }

    return Iterator {
      std::begin(container),
      std::move(maybeNextToBegin)
    };
  }

  Iterator end() const {
    const auto& container = ContainerBinding::value;

    auto maybePriorToEnd = std::end(container);
    if(maybePriorToEnd != std::begin(container)) {
      --maybePriorToEnd;
    }

    return Iterator {
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
} // namespace Scine

#endif
