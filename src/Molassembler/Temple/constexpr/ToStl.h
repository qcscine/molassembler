/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Create STL analogues of constexpr containers
 *
 * Offers quick functions to create STL analog data structures from the ones
 * defined in this library.
 */

#ifndef INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_TO_STL_H
#define INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_TO_STL_H

#include "Molassembler/Temple/constexpr/Array.h"
#include "Molassembler/Temple/constexpr/DynamicMap.h"

#include <vector>
#include <map>
#include <set>

namespace Scine {
namespace Molassembler {
namespace Temple {

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

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
