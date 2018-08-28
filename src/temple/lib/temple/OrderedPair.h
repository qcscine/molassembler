#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ORDERED_PAIR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ORDERED_PAIR_H

#include <tuple>
#include <utility>

/*!@file
 *
 * Contains a class imitating std::pair whose member types are homogeneous and
 * ordered.
 */

/* TODO
 * - Iterator testing
 */

namespace temple {

template<typename T>
struct OrderedPair {
  // Standard guarantees first and second are laid out sucessively in memory
  T first;
  T second;

  OrderedPair() = default;

  constexpr OrderedPair(T a, T b) : first {std::move(a)}, second {std::move(b)} {
    if(b < a) {
      std::swap(first, second);
    }
  }

  constexpr T front() const {
    return first;
  }

  constexpr T& front() {
    return first;
  }

  constexpr T back() const {
    return second;
  }

  constexpr T& back() {
    return second;
  }

  constexpr T& at(unsigned i) {
    if(i > 1) {
      throw std::out_of_range("Invalid access to pair");
    }

    if(i == 0) {
      return first;
    }

    return second;
  }

  constexpr T at(unsigned i) const {
    if(i > 1) {
      throw std::out_of_range("Invalid access to pair");
    }

    if(i == 0) {
      return first;
    }

    return second;
  }

  constexpr T& operator [] (unsigned i) {
    return at(i);
  }

  constexpr T operator [] (unsigned i) const {
    return at(i);
  }

  class Iterator {
  public:
    // TODO iterator typedefs
    /*using iterator_category = std::bidirectional_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = T&;*/

  private:
    T* _ptr;

  public:
    Iterator() = default;
    Iterator(T* ptr) : _ptr {ptr} {}

    Iterator& operator ++ () {
      ++_ptr;
      return *this;
    }

    Iterator operator ++ (int) {
      Iterator retval = *this;
      ++(*this);
      return retval;
    }

    Iterator& operator -- () {
      --_ptr;
      return *this;
    }

    Iterator operator -- (int) {
      Iterator retval = *this;
      --(*this);
      return retval;
    }

    T& operator * () const {
      return *_ptr;
    }

    std::ptrdiff_t operator - (const Iterator& other) const {
      return other._ptr - _ptr;
    }

    bool operator == (const Iterator& other) const {
      return _ptr == other._ptr;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }
  };

  Iterator begin() {
    return Iterator {&first};
  }

  Iterator end() {
    // Make a past-the-end iterator
    Iterator iter {&second};
    ++iter;
    return iter;
  }

  class ConstIterator {
  public:
    // TODO iterator typedefs

  private:
    T* const _ptr;

  public:
    ConstIterator() = default;
    ConstIterator(T* ptr) : _ptr {ptr} {}

    ConstIterator& operator ++ () {
      ++_ptr;
      return *this;
    }

    ConstIterator operator ++ (int) {
      ConstIterator retval = *this;
      ++(*this);
      return retval;
    }

    ConstIterator& operator -- () {
      --_ptr;
      return *this;
    }

    ConstIterator operator -- (int) {
      ConstIterator retval = *this;
      --(*this);
      return retval;
    }

    const T& operator * () const {
      return *_ptr;
    }

    std::ptrdiff_t operator - (const Iterator& other) const {
      return other._ptr - _ptr;
    }

    bool operator == (const ConstIterator& other) const {
      return _ptr == other._ptr;
    }

    bool operator != (const ConstIterator& other) const {
      return !(*this == other);
    }
  };

  ConstIterator begin() const {
    return ConstIterator {&first};
  }

  ConstIterator end() const {
    // Make a past-the-end iterator
    ConstIterator iter {&second};
    ++iter;
    return iter;
  }

  // C++17 spaceship operator
  constexpr bool operator < (const OrderedPair& other) const {
    return std::tie(first, second) < std::tie(other.first, other.second);
  }

  constexpr bool operator > (const OrderedPair& other) const {
    return std::tie(first, second) > std::tie(other.first, other.second);
  }

  constexpr bool operator == (const OrderedPair& other) const {
    return std::tie(first, second) == std::tie(other.first, other.second);
  }

  constexpr bool operator != (const OrderedPair& other) const {
    return std::tie(first, second) != std::tie(other.first, other.second);
  }

  template<typename UnaryFunction>
  auto map(UnaryFunction&& mapFunction) const {
    return std::make_pair(
      mapFunction(first),
      mapFunction(second)
    );
  }
};

} // namespace temple

#endif
