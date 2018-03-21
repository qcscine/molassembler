#ifndef INCLUDE_CONSTEXPR_MAGIC_TO_STL_H
#define INCLUDE_CONSTEXPR_MAGIC_TO_STL_H

#include "Array.h"
#include "DynamicMap.h"

#include <vector>
#include <map>
#include <set>

/*! @file
 *
 * Offers quick functions to create STL analog data structures from the ones
 * defined in this library.
 */

namespace temple {

//! Converts an Array into a std::array
template<typename T, size_t size>
std::array<T, size> toSTL(const Array<T, size>& array) {
  return array.getArray();
}

//! Converts a DynamicArray into a std::vector
template<typename T, size_t size>
std::vector<T> toSTL(const DynamicArray<T, size>& dynamicArray) {
  std::vector<T> data;
  data.reserve(dynamicArray.size());

  for(unsigned i = 0; i < dynamicArray.size(); ++i) {
    data.push_back(dynamicArray.at(i));
  }

  return data;
}

//! Converts a DynamicSet into a std::set
template<typename T, size_t size>
std::set<T> toSTL(const DynamicSet<T, size>& dynamicSet) {
  std::set<T> returnSet;

  for(const auto& element : dynamicSet) {
    returnSet.insert(element);
  }

  return returnSet;
}

/*template<typename T, typename U, size_t size>
std::map<T, U> toSTL(const DynamicMap<T, U, size>& dynamicMap) {
  std::map<T, U> returnMap;

  for(
}*/

} // namespace temple

#endif
