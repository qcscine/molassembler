/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides enumeration for easy range-for use of unindexed containers
 *
 * Provides a range-for compatible struct exposing begin and end forward
 * iterators that have a pair as their value_type, the first of which is the
 * current index, the second of which is the current object in the container.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ENUMERATE_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ENUMERATE_ADAPTOR_H

#include "molassembler/Temple/Preprocessor.h"
#include "molassembler/Temple/Binding.h"

#include <memory>

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Adaptors {
namespace Detail {

template<class Container>
class Enumerator : public Binding<Container> {
public:
  using ContainerBinding = Binding<Container>;

  // Get the bare Container containing type
  using T = decltype(
    *std::begin(
      std::declval<const Container>()
    )
  );

  using ContainerIterator = decltype(
    std::begin(
      std::declval<const Container>()
    )
  );

  struct EnumerationStruct {
    const unsigned index;
    T value;

    EnumerationStruct(
      const unsigned passIndex,
      const T& passValue
    ) : index(passIndex),
        value(passValue)
    {}
  };

  using ContainerBinding::ContainerBinding;

  class iterator {
  private:
    ContainerIterator it_;
    unsigned index_;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = EnumerationStruct;
    using difference_type = unsigned;
    using pointer = const EnumerationStruct*;
    using reference = const EnumerationStruct;

    explicit iterator(
      ContainerIterator it,
      unsigned index
    ) : it_(std::move(it)),
        index_(index)
    {}

    iterator& operator ++ () {
      ++it_;
      ++index_;
      return *this;
    }

    iterator operator++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    bool operator == (iterator other) const {
      return it_ == other.it_;
    }

    bool operator != (iterator other) const {
      return !(*this == other);
    }

    reference operator * () const {
      return EnumerationStruct {
        index_,
        *it_
      };
    }
  };

  iterator begin() const {
    return iterator(
      std::begin(ContainerBinding::value),
      0
    );
  }

  iterator end() const {
    return iterator(
      std::end(ContainerBinding::value),
      ContainerBinding::value.size()
    );
  }
};

} // namespace Detail

/*! Returns an EnerateTemporary for use with range-for expressions that
 * generates a struct with members index and value for every contained element.
 * Requires that the container implements begin(), end() and size() members.
 * Should involve minimal copying, mostly uses references, though no space or
 * time complexity guarantees are given.
 *
 * This may seem like overkill, particularly when most containers' elements can
 * be accessed with operator [], where customary index-based loops are more
 * adequate. However, perhaps some custom containers do not implement operator
 * [], and for these, use this.
 */
template<class Container>
Detail::Enumerator<Container> enumerate(Container&& container) {
  return Detail::Enumerator<Container>(
    std::forward<Container>(container)
  );
}

} // namespace Adaptors
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
