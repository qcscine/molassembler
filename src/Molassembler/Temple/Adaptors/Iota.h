/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides an integer range generator.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_IOTA_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_IOTA_ADAPTOR_H

#include "Molassembler/Temple/Traits.h"

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Adaptors {
namespace Detail {

template<typename T>
class Iota {
public:
//!@name Constructors
//!@{
  Iota(T lower, T upper) : lower_(lower), upper_(upper) {
    assert(lower_ < upper_);
  }

  explicit Iota(T upper) : Iota(T {0}, upper) {}
//!@}

  std::size_t size() const {
    // This is fine since we know that upper_ > lower_ (see ctor)
    return upper_ - lower_;
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
    return iterator {lower_};
  }

  iterator end() const {
    return iterator {upper_};
  }
//!@}

private:
//!@name State
//!@{
  T lower_, upper_;
//!@}
};

} // namespace Detail

template<typename T>
auto range(T lower, T upper) {
  return Detail::Iota<T>(lower, upper);
}

template<typename T>
auto range(T upper) {
  return Detail::Iota<T>(upper);
}

} // namespace Adaptors
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
