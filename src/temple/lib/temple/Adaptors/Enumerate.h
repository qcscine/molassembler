#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ENUMERATE_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ENUMERATE_ADAPTOR_H

#include "temple/Preprocessor.h"

#include <memory>

/*! @file
 *
 * @brief Provides enumeration for easy range-for use of unindexed containers
 *
 * Provides a range-for compatible struct exposing begin and end forward
 * iterators that have a pair as their value_type, the first of which is the
 * current index, the second of which is the current object in the container.
 */

namespace temple {

namespace adaptors {

namespace detail {

template<class Container>
class Enumerator {
public:
  // See tricks documentation
  using BoundContainer = std::conditional_t<
    std::is_rvalue_reference<Container&&>::value,
    std::decay_t<Container>,
    const Container&
  >;

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

private:
  BoundContainer _container;

public:
  Enumerator(Container&& container)
    : _container(std::forward<Container>(container)) {}

  class iterator {
  private:
    ContainerIterator _it;
    unsigned _index;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = EnumerationStruct;
    using difference_type = unsigned;
    using pointer = const EnumerationStruct*;
    using reference = const EnumerationStruct;

    explicit iterator(
      ContainerIterator it,
      unsigned index
    ) : _it(std::move(it)),
        _index(index)
    {}

    iterator& operator ++ () {
      ++_it;
      ++_index;
      return *this;
    }

    iterator operator++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    bool operator == (iterator other) const {
      return _it == other._it;
    }

    bool operator != (iterator other) const {
      return !(*this == other);
    }

    reference operator * () const {
      return EnumerationStruct {
        _index,
        *_it
      };
    }
  };

  iterator begin() const {
    return iterator(
      std::begin(_container),
      0
    );
  }

  iterator end() const {
    return iterator(
      std::end(_container),
      _container.size()
    );
  }
};

} // namespace detail

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
detail::Enumerator<Container> enumerate(Container&& container) {
  return detail::Enumerator<Container>(
    std::forward<Container>(container)
  );
}

} // namespace adaptors

} // namespace temple

#endif
