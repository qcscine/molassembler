#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_IOTA_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_IOTA_ADAPTOR_H

#include "temple/Traits.h"

namespace temple {

namespace adaptors {

namespace detail {

template<typename T>
class Iota {
public:
//!@name Constructors
//!@{
  Iota(T lower, T upper) : _lower(lower), _upper(upper) {
    static_assert(
      std::is_integral<T>::value,
      "Iota template argument must be an integral value"
    );

    assert(_lower < _upper);
  }

  explicit Iota(T upper) : Iota(T {0}, upper) {}
//!@}

  std::size_t size() const {
    // This is fine since we know that _upper > _lower (see ctor)
    return _upper - _lower;
  }

//!@name Iterators
//!@{
  struct iterator {
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T;
    using reference = T;

    T value;

    explicit iterator() = default;
    explicit iterator(T initValue) : value(initValue) {}

    // Prefix increment
    iterator& operator ++ () {
      ++value;
      return *this;
    }

    // Postfix increment
    iterator operator ++ (int) {
      iterator prior = *this;
      ++(*this);
      return prior;
    }

    iterator operator + (const T increment) const {
      return iterator {
        value + increment
      };
    }

    iterator& operator += (const T increment) {
      value += increment;
      return *this;
    }

    iterator operator - (const T decrement) const {
      return iterator {
        value - decrement
      };
    }

    iterator& operator -= (const T decrement) {
      value -= decrement;
      return *this;
    }

    iterator& operator -- () {
      --value;
      return *this;
    }

    iterator operator -- (int) {
      iterator prior = *this;
      --(*this);
      return prior;
    }

    std::ptrdiff_t operator - (const iterator& other) const {
      return (
        static_cast<std::ptrdiff_t>(value)
        - static_cast<std::ptrdiff_t>(other.value)
      );
    }

    T operator * () const {
      return value;
    }

    bool operator == (const iterator& other) const {
      return value == other.value;
    }

    bool operator != (const iterator& other) const {
      return !(*this == other);
    }
  };

  iterator begin() const {
    return iterator {_lower};
  }

  iterator end() const {
    return iterator {_upper};
  }
//!@}

private:
//!@name State
//!@{
  T _lower, _upper;
//!@}
};

} // namespace detail

template<typename IntegerType>
auto range(IntegerType lower, IntegerType upper) {
  static_assert(
    std::is_integral<IntegerType>::value,
    "Range must be called with an integer value"
  );

  return detail::Iota<IntegerType>(lower, upper);
}

template<typename IntegerType>
auto range(IntegerType upper) {
  static_assert(
    std::is_integral<IntegerType>::value,
    "Range must be called with an integer value"
  );

  return detail::Iota<IntegerType>(upper);
}

} // namespace adaptors

} // namespace temple


#endif
