#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_UTIL_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_UTIL_H

#include <set>

#include "boost/optional.hpp"

/*!@file
 *
 * @brief Some few still-used utility functions.
 */

/* TODO
 * - phase these out from the few places where they are in use in this library
 */

namespace stereopermutation {

namespace Util {

template<typename Function, typename T1, typename T2>
decltype(auto) minMaxAdapt(
  Function function,
  T1 a,
  T2 b
) {
  return function(
    std::min(a, b),
    std::max(a, b)
  );
}

template<typename Function, typename T>
decltype(auto) sortBinaryArgs(
  Function&& function,
  T&& a,
  T&& b
) {
  if(a < b) {
    return function(std::forward<T>(a), std::forward<T>(b));
  }

  return function(std::forward<T>(b), std::forward<T>(a));
}

} // namespace Util

} // namespace stereopermutation

#endif
