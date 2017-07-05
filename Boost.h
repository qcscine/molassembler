#ifndef INCLUDE_TEMPLATE_MAGIC_BOOST_ENHANCEMENTS_H
#define INCLUDE_TEMPLATE_MAGIC_BOOST_ENHANCEMENTS_H

#include <boost/optional.hpp>

namespace TemplateMagic {

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
boost::optional<bool> componentSmaller(
  const T& a,
  const T& b
) {
  if(a < b) return true;
  else if(b < a) return false;
  else return boost::none;
}

} // namespace TemplateMagic

#endif
