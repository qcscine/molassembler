#ifndef INCLUDE_CONSTEXPR_MAGIC_BOOST_H
#define INCLUDE_CONSTEXPR_MAGIC_BOOST_H

#include "Optional.h"

namespace constable {

/*! 
 * To help with creating consistent logical operators for multi-component
 * structs / classes. E.g.::
 *   
 *   struct Foo {
 *     unsigned a, b, c;
 *
 *     bool operator < (const Foo& other) {
 *       return componentSmaller(a, other.a).value_or( // sort by a first
 *         componentSmaller(b, other.b).value_or( // then by b
 *           componentSmaller(c, other.c).value_or( // then by c
 *             false
 *           )
 *         )
 *       );
 *     }
 *
 *     // the alternative can quickly become a lengthy chain and does not 
 *     // terminate as quickly if x > other.x
 *     bool operator < (const Foo& other) {
 *       return (
 *         a < other.a
 *         || (
 *           a == other.a
 *           && b < other.b
 *         ) || (
 *           a == other.a
 *           && b == other.b
 *           && c < other.c
 *         )
 *       );
 *     }
 *   };
 *
 */
template<typename T>
constexpr Optional<bool> componentSmaller(
  const T& a,
  const T& b
) {
  if(a < b) return true;
  else if(b < a) return false;
  else return {};
}

} // namespace constable

#endif
