/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Helper functions for dealing with optional types
 *
 * Implements optional-returning function composition syntactic sugar to avoid
 * repetitive patterns when dealing with lots of optionals.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_OPTIONALS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_OPTIONALS_H

#include <boost/optional.hpp>
#include "temple/Traits.h"

namespace temple {

namespace optionals {

// UnaryFunction: T -> U
template<
  template<typename> class OptionalType,
  typename T,
  class UnaryFunction
>
auto map(const OptionalType<T>& optional, UnaryFunction&& function) {
  using U = decltype(function(std::declval<T>()));

  if(optional) {
    return OptionalType<U> {function(optional.value())};
  }

  return OptionalType<U> {};
}

// UnaryFunction: T -> Optional<U>
template<
  template<typename> class OptionalType,
  typename T,
  class UnaryFunction
>
auto flatMap(const OptionalType<T>& optional, UnaryFunction&& function) {
  using OptionalU = decltype(function(std::declval<T>()));

  if(optional) {
    return OptionalU {function(optional.value())};
  }

  return OptionalU {};
}

} // namespace optionals

} // namespace temple

#endif
