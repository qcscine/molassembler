/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides functional-style transformation of container elements
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_TRANSFORM_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_TRANSFORM_ADAPTOR_H

#include "Molassembler/Temple/ContainerTraits.h"
#include "Molassembler/Temple/Invoke.h"
#include "Molassembler/Temple/Binding.h"

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Adaptors {
namespace Detail {

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
    Temple::invoke(
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
    Traits::hasSize<Container>::value,
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
      return Temple::invoke(basePtr_->function, *iter_);
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

} // namespace Detail

template<class Container, typename UnaryFunction>
auto transform(
  Container&& container,
  UnaryFunction&& function
) {
  return Detail::Transformer<Container, UnaryFunction>(
    std::forward<Container>(container),
    std::forward<UnaryFunction>(function)
  );
}

} // namespace Adaptors
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
