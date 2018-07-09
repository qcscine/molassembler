#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ORDERED_PAIR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ORDERED_PAIR_H

#include <tuple>

namespace temple {

template<typename T>
struct OrderedPair {
  // Standard guarantees first and second are laid out sucessively in memory
  T first;
  T second;

  OrderedPair() = default;
  OrderedPair(T a, T b) : first {std::min(a, b)}, second {std::max(a, b)} {}

  T front() const {
    return first;
  }

  T& front() {
    return first;
  }

  T back() const {
    return second;
  }

  T& back() {
    return second;
  }

  T& at(unsigned i) {
    if(i > 1) {
      throw std::out_of_range("Invalid access to pair");
    }

    if(i == 0) {
      return first;
    }

    return second;
  }

  T at(unsigned i) const {
    if(i > 1) {
      throw std::out_of_range("Invalid access to pair");
    }

    if(i == 0) {
      return first;
    }

    return second;
  }

  T& operator [] (unsigned i) {
    return at(i);
  }

  T operator [] (unsigned i) const {
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
};

} // namespace temple

#endif
