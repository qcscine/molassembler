#ifndef INCLUDE_CONSTEXPR_MAGIC_PAIR_H
#define INCLUDE_CONSTEXPR_MAGIC_PAIR_H

#include <utility>

/*! @file
 *
 * Implements a constexpr container much like std::pair.
 */

namespace ConstexprMagic {

/*! 
 * Heterogeneous pair type.
 *
 * Requires that both types T and U are default-constructible.
 */
template<typename T, typename U>
struct Pair {
  T first;
  U second;

  constexpr Pair() : first(T {}), second (U {}) {}

  /* Value initialization */
  constexpr Pair(const T& passFirst, const U& passSecond) 
    : first(passFirst), second(passSecond)
  {}

  constexpr Pair(T&& passFirst, U&& passSecond) 
    : first(passFirst), second(passSecond)
  {}

  /* Copy initialization */
  constexpr Pair(const Pair& other) 
    : first(other.first), second(other.second) 
  {}

  constexpr Pair(Pair&& other) 
    : first(other.first), second(other.second) 
  {}

  /* Assignment */
  constexpr Pair& operator = (const Pair& other) {
    first = other.first;
    second = other.second;

    return *this;
  }

  constexpr Pair& operator = (Pair&& other) {
    first = std::move(other.first);
    second = std::move(other.second);

    return *this;
  }

  /* Comparison operators */
  constexpr bool operator < (const Pair& other) const {
    if(first > other.first) {
      return false;
    } 
    
    if(first == other.first) {
      return second < other.second;
    }

    return true;
  }

  constexpr bool operator > (const Pair& other) const {
    return (other < *this);
  }

  constexpr bool operator == (const Pair& other) const {
    return (
      first == other.first
      && second == other.second
    );
  }

  constexpr bool operator != (const Pair& other) const {
    return !(*this == other);
  }
};

} // namespace ConstexprMagic

#endif
