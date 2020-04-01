/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides functional-style transformation of container elements
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_TRANSFORM_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_TRANSFORM_ADAPTOR_H

#include "molassembler/Temple/ContainerTraits.h"
#include "molassembler/Temple/Invoke.h"
#include "molassembler/Temple/Binding.h"

namespace Scine {
namespace temple {
namespace adaptors {
namespace detail {

template<class Container, typename UnaryFunction>
struct Transformer {
//!@name Types
//!@{
  using BoundContainer = typename Binding<Container>::type;

  using ContainerValueType = decltype(
    *std::begin(
      std::declval<const Container>()
    )
  );

  using ReturnType = decltype(
    invoke(
      std::declval<UnaryFunction>(),
      std::declval<ContainerValueType>()
    )
  );

  using ContainerIteratorType = decltype(std::begin(std::declval<const Container>()));
//!@}

//!@name State
//!@{
  BoundContainer container;
  UnaryFunction function;
//!@}

//!@name Constructor
//!@{
  Transformer(
    Container&& passContainer,
    UnaryFunction&& passFunction
  ) : container(std::forward<Container>(passContainer)),
      function(passFunction)
  {}
//!@}

//!@name Information
//!@{
  std::enable_if_t<
    traits::hasSize<Container>::value,
    std::size_t
  > size() const {
    return container.size();
  }
//!@}

//!@name Iterators
//!@{
  class iterator {
  private:
    const Transformer* basePtr_;
    ContainerIteratorType iter_;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ReturnType;
    using difference_type = int;
    using pointer = const ReturnType*;
    using reference = const ReturnType&;

    iterator() = default;
    iterator(
      const Transformer& base,
      ContainerIteratorType&& iter
    ) : basePtr_(&base),
        iter_ {iter}
    {}

    iterator& operator ++ () {
      ++iter_;
      return *this;
    }

    iterator operator ++ (int) {
      iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const iterator& other) const {
      return iter_ == other.iter_;
    }

    bool operator != (const iterator& other) const {
      return !(*this == other);
    }

    ReturnType operator * () const {
      return invoke(basePtr_->function, *iter_);
    }
  };

  iterator begin() const {
    return {
      *this,
      std::begin(container)
    };
  }

  iterator end() const {
    return {
      *this,
      std::end(container)
    };
  }
//!@}
};

} // namespace detail

template<class Container, typename UnaryFunction>
auto transform(
  Container&& container,
  UnaryFunction&& function
) {
  return detail::Transformer<Container, UnaryFunction>(
    std::forward<Container>(container),
    std::forward<UnaryFunction>(function)
  );
}

} // namespace adaptors
} // namespace temple
} // namespace Scine

#endif
