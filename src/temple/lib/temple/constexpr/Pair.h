#ifndef INLCUDE_MOLASSEMBLER_TEMPLE_PAIR_H
#define INLCUDE_MOLASSEMBLER_TEMPLE_PAIR_H

#include "temple/Preprocessor.h"

#include <utility>

/*! @file
 *
 * @brief std::pair-like container
 *
 * Implements a constexpr container much like std::pair. Yes, this really is
 * necessary since C++14's pair isn't constexpr everywhere!
 */

namespace temple {

/*!
 * Heterogeneous pair type.
 *
 * Requires that both types T and U are default-constructible.
 */
template<typename T, typename U>
struct Pair {
//!@name State
//!@{
  T first;
  U second;
//!@}

//!@name Constructors
//!@{
  constexpr Pair() : first(T {}), second (U {}) {}

  /* Value initialization */
  constexpr Pair(const T& passFirst, const U& passSecond)
    : first(passFirst), second(passSecond)
  {}

  constexpr Pair(T&& passFirst, U&& passSecond)
    : first(passFirst), second(passSecond)
  {}
//!@}

//!@name Special member functions
  constexpr Pair(const Pair& other)
    : first(other.first), second(other.second) {}
  constexpr Pair(Pair&& other) noexcept
    : first(other.first), second(other.second) {}
  constexpr Pair& operator = (const Pair& other) {
    first = other.first;
    second = other.second;

    return *this;
  }
  constexpr Pair& operator = (Pair&& other) noexcept {
    first = std::move(other.first);
    second = std::move(other.second);

    return *this;
  }
  ~Pair() = default;
//!@}

//!@name Operators
//!@{
  //! Lexicographical comparison
  PURITY_WEAK constexpr bool operator < (const Pair& other) const {
    if(first > other.first) {
      return false;
    }

    if(first == other.first) {
      return second < other.second;
    }

    return true;
  }

  //! Lexicographical comparison
  PURITY_WEAK constexpr bool operator > (const Pair& other) const {
    return (other < *this);
  }

  //! Lexicographical comparison
  PURITY_WEAK constexpr bool operator == (const Pair& other) const {
    return (
      first == other.first
      && second == other.second
    );
  }

  //! Lexicographical comparison
  PURITY_WEAK constexpr bool operator != (const Pair& other) const {
    return !(*this == other);
  }
//!@}
};

} // namespace temple

#endif
