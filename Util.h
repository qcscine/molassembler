#ifndef LIB_UNIQUE_ASSIGNMENTS_UTIL_H
#define LIB_UNIQUE_ASSIGNMENTS_UTIL_H

#include <sstream>
#include <set>

#include "boost/optional.hpp"

namespace UniqueAssignments {

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

template<typename T>
std::string toString(const T& container) {
  std::stringstream sstream;
  sstream << "{";
  unsigned nItems = container.size();
  for(const auto& item : container) {
    sstream << item;
    if(--nItems != 0) sstream << ", ";
  }
  sstream << "}";

  return sstream.str();
}

} // eo namespace Util

} // eo namespace UniqueAssignments

#endif
