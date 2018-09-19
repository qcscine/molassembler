// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_TRANSFORM_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_TRANSFORM_ADAPTOR_H

#include "temple/ContainerTraits.h"
#include "temple/Invoke.h"

/*!@file
 *
 * @brief Provides functional-style transformation of container elements
 */

namespace temple {

namespace adaptors {

namespace detail {

template<class Container, typename UnaryFunction>
struct Transformer {
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
    const Transformer* _basePtr;
    ContainerIteratorType _iter;

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
    ) : _basePtr(&base),
        _iter {iter}
    {}

    iterator& operator ++ () {
      ++_iter;
      return *this;
    }

    iterator operator ++ (int) {
      iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const iterator& other) const {
      return _iter == other._iter;
    }

    bool operator != (const iterator& other) const {
      return !(*this == other);
    }

    ReturnType operator * () const {
      return invoke(_basePtr->function, *_iter);
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

#endif
