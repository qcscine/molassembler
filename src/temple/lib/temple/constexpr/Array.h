#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ARRAY_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ARRAY_H

#include "temple/constexpr/Containers.h"

#include <cstddef>
#include <type_traits>
#include <utility>

/*! @file
 *
 * @brief std::array-like class
 *
 * Constexpr fixed-size array to replace std::array in C++14. This class is
 * largely unneeded in C++17 since many std::array members are then marked
 * constexpr.
 */

namespace temple {

template<typename T, size_t nItems>
class Array {
public:
//!@name Types
//!@{
  using value_type = T;
//!@}

//!@name Constructors
//!@{
  /*! Delegate constructor using another Array and an index_sequence to directly
   * form the array mem-initializer with a parameter pack expansion
   */
  template<size_t ... Inds>
  constexpr Array(
    const Array& other,
    std::index_sequence<Inds...> /* inds */
  )
    :_items {other[Inds]...}
  {}

  /* Constructing from another array is tricky since we're technically not
   * allowed to edit _items in-class, so we delegate to the previous constructor
   * and directly form the mem-initializer
   */
  constexpr Array(const Array& other)
    : Array(other, std::make_index_sequence<nItems>{})
  {}

  //! Delegate std::array ctor, using same trick as copy ctor
  template<size_t ... Inds>
  constexpr Array(
    const std::array<T, nItems>& other,
    std::index_sequence<Inds...> /* inds */
  ) :_items {other[Inds]...} {}

  //! Construct from std::array using same trick as copy ctor
  constexpr Array(const std::array<T, nItems>& other)
    : Array(other, std::make_index_sequence<nItems>{})
  {}

  //! Parameter pack constructor, will work as long as the arguments are castable
  template<typename ...Args>
  constexpr Array(Args... args)
    : _items {static_cast<T>(args)...}
  {}
//!@}

  constexpr T& operator[] (const std::size_t index) noexcept {
    return _items[index];
  }

  constexpr T& at(const std::size_t index) noexcept {
    return _items[index];
  }

  constexpr const T& operator[] (const std::size_t index) const noexcept {
    return _items[index];
  }

  constexpr const T& at(const std::size_t index) const noexcept {
    return _items[index];
  }

  constexpr const T& front() const noexcept {
    return _items[0];
  }

  constexpr const T& back() const noexcept {
    return _items[nItems - 1];
  }

  constexpr size_t size() const noexcept {
    return nItems;
  }

  constexpr bool operator == (const Array& other) const noexcept {
    for(std::size_t i = 0; i < nItems; ++i) {
      if(_items[i] != other._items[i]) {
        return false;
      }
    }

    return true;
  }

  constexpr bool operator != (const Array& other) const noexcept {
    for(std::size_t i = 0; i < nItems; ++i) {
      if(_items[i] != other._items[i]) {
        return true;
      }
    }

    return false;
  }

  constexpr bool operator < (const Array& other) const noexcept {
    for(std::size_t i = 0; i < nItems; ++i) {
      if(_items[i] < other._items[i]) {
        return true;
      }

      if(_items[i] > other._items[i]) {
        return false;
      }
    }

    return false;
  }

  constexpr bool operator > (const Array& other) const noexcept {
    return (other < *this);
  }

  // Begin and end iterators
  class iterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = int;
    using pointer = const T*;
    using reference = T&;

    constexpr iterator(
      Array& instance,
      std::size_t initPosition
    ) noexcept
      : _baseRef(instance),
        _position(initPosition)
    {}

    constexpr iterator() = delete;
    constexpr iterator(const iterator& other) noexcept
      : _baseRef(other._baseRef),
        _position(other._position)
    {}
    constexpr iterator(iterator&& other) noexcept = default;
    constexpr iterator& operator = (const iterator& other) noexcept {
      _baseRef = other._baseRef;
      _position = other._position;

      return *this;
    }
    constexpr iterator& operator = (iterator&& other) noexcept = default;
    ~iterator() = default;

    constexpr iterator& operator ++ () noexcept {
      _position += 1;
      return *this;
    }

    constexpr iterator operator ++ (int) noexcept {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr iterator& operator --() noexcept {
      _position -= 1;
      return *this;
    }

    constexpr iterator operator -- (int) noexcept {
      iterator retval = *this;
      --(*this);
      return retval;
    }

    constexpr iterator operator + (const int increment) noexcept {
      iterator retval = *this;
      retval += increment;
      return retval;
    }

    constexpr iterator operator - (const int increment) noexcept {
      iterator retval = *this;
      retval -= increment;
      return retval;
    }

    constexpr iterator& operator += (const int increment) noexcept {
      _position += increment;
      return *this;
    }

    constexpr iterator& operator -= (const int increment) noexcept {
      _position -= increment;
      return *this;
    }

    PURITY_WEAK constexpr int operator - (const iterator& other) const noexcept {
      return (
        static_cast<int>(_position)
        - static_cast<int>(other._position)
      );
    }

    PURITY_WEAK constexpr bool operator == (const iterator& other) const noexcept {
      return (
        &_baseRef == &other._baseRef
        && _position == other._position
      );
    }

    PURITY_WEAK constexpr bool operator != (const iterator& other) const noexcept {
      return !(
        *this == other
      );
    }

    PURITY_WEAK constexpr reference operator * () const noexcept {
      return _baseRef[_position];
    }

  private:
    Array& _baseRef;
    std::size_t _position;
  };

  PURITY_WEAK constexpr iterator begin() noexcept {
    return iterator(*this, 0);
  }

  PURITY_WEAK constexpr iterator end() noexcept {
    return iterator(*this, nItems);
  }

  class const_iterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = int;
    using pointer = const T*;
    using reference = const T&;

    constexpr explicit const_iterator(
      const Array& instance,
      std::size_t initPosition
    ) noexcept
      : _baseRef(instance),
        _position(initPosition)
    {}

    constexpr const_iterator() = delete;
    constexpr const_iterator(const const_iterator& other) noexcept
      : _baseRef(other._baseRef),
        _position(other._position)
    {}
    constexpr const_iterator(const_iterator&& other) noexcept = default;
    constexpr const_iterator& operator = (const const_iterator& other) {
      if(_baseRef != other._baseRef) {
        throw "Trying to assign constIterator to other base Array!";
      }

      _position = other._position;

      return *this;
    }
    constexpr const_iterator& operator = (const_iterator&& other) noexcept = default;
    ~const_iterator() = default;

    constexpr const_iterator& operator ++ () noexcept {
      _position += 1;
      return *this;
    }

    constexpr const_iterator operator ++ (int) noexcept {
      const_iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr const_iterator& operator --() noexcept {
      _position -= 1;
      return *this;
    }

    constexpr const_iterator operator -- (int) noexcept {
      const_iterator retval = *this;
      --(*this);
      return retval;
    }

    constexpr const_iterator operator + (const int increment) noexcept {
      const_iterator retval = *this;
      retval += increment;
      return retval;
    }

    constexpr const_iterator operator - (const int increment) noexcept {
      const_iterator retval = *this;
      retval -= increment;
      return retval;
    }

    constexpr const_iterator& operator += (const int increment) noexcept {
      _position += increment;
      return *this;
    }

    constexpr const_iterator& operator -= (const int increment) noexcept {
      _position -= increment;
      return *this;
    }

    constexpr int operator - (const const_iterator& other) const noexcept {
      return (
        static_cast<int>(_position)
        - static_cast<int>(other._position)
      );
    }

    constexpr bool operator == (const const_iterator& other) const noexcept {
      return (
        &_baseRef == &other._baseRef
        && _position == other._position
      );
    }

    constexpr bool operator != (const const_iterator& other) const noexcept {
      return !(
        *this == other
      );
    }

    constexpr reference operator * () const noexcept {
      return _baseRef[_position];
    }

  private:
    const Array& _baseRef;
    std::size_t _position;
  };

  constexpr const_iterator begin() const noexcept {
    return const_iterator(*this, 0);
  }

  constexpr const_iterator end() const noexcept {
    return const_iterator(*this, nItems);
  }

  //! Implicit conversion operator to a std::array
  constexpr operator std::array<T, nItems> () const {
    return _makeArray(std::make_index_sequence<nItems>{});
  }

  // Explicit conversion to a std::array
  constexpr std::array<T, nItems> getArray() const {
    return _makeArray(std::make_index_sequence<nItems>{});
  }

private:
  T _items[nItems];

  template<size_t ... Inds>
  std::array<T, nItems> _makeArray(std::index_sequence<Inds...> /* inds */) {
    return {{
      _items[Inds]...
    }};
  }
};

// This way, makeArray must be called with at least one argument
template<typename T, typename... Tail>
constexpr auto makeArray(T head, Tail... tail) -> Array<T, 1 + sizeof...(Tail)> {
  return { head, tail ... };
}

} // namespace temple

#endif
