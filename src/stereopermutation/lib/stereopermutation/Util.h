#ifndef LIB_UNIQUE_ASSIGNMENTS_UTIL_H
#define LIB_UNIQUE_ASSIGNMENTS_UTIL_H

#include <set>

#include "boost/optional.hpp"

namespace stereopermutation {

namespace Util {

template<typename Comparable>
boost::optional<bool> compareSmaller(
  const Comparable& a,
  const Comparable& b
) {
  if(a < b) return true;
  else if(b < a) return false;
  else return {};
}

template<typename T>
std::set<T> removeElementFromSet(
  std::set<T> set,
  T toRemove
) {
  set.erase(toRemove);
  return set;
}

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
  } else {
    return function(std::forward<T>(b), std::forward<T>(a));
  }
}

} // eo namespace Util

} // eo namespace stereopermutation

#endif
