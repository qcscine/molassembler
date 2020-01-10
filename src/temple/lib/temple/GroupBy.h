/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Group container elements together by various methods
 *
 * Provides functionality to group elements of a container together by various
 * methods.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_GROUP_BY_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_GROUP_BY_H

#include "temple/Traits.h"

#include <vector>
#include <map>

namespace Scine {
namespace temple {

/*!
 * Split a container's values by a mapping function whose return value elements
 * are compared with. Matching mapped values are grouped and returned in a
 * ragged 2D vector.
 */
template<class Container, class UnaryFunction>
std::vector<
  std::vector<
    traits::getValueType<Container>
  >
> groupByMapping(
  const Container& container,
  UnaryFunction&& function
) {
  using T = traits::getValueType<Container>;
  using R = traits::functionReturnType<UnaryFunction, T>;

  std::vector<
    std::vector<T>
  > groups;

  std::map<R, unsigned> indexMap;

  for(auto iter = std::begin(container); iter != std::end(container); ++iter) {
    auto ret = function(*iter);
    if(indexMap.count(ret) == 0) {
      indexMap[ret] = groups.size();

      groups.emplace_back(
        std::vector<T> {*iter}
      );
    } else {
      groups.at(
        indexMap.at(ret)
      ).push_back(*iter);
    }
  }

  return groups;
}

/*!
 * Split a container's values by a binary comparison function. Returns a ragged
 * 2D vector.
 *
 * @note Requires that the equality comparison is symmetric and transitive!
 */
template<class Container, class BinaryFunction>
std::vector<
  std::vector<
    traits::getValueType<Container>
  >
> groupByEquality(
  const Container& container,
  BinaryFunction&& compareEqual
) {
  using T = traits::getValueType<Container>;

  std::vector<
    std::vector<T>
  > groups;

  for(auto iter = std::begin(container); iter != std::end(container); ++iter) {
    bool foundEqual = false;
    for(auto& group : groups) {

      if(compareEqual(*iter, *std::begin(group))) {
        group.push_back(*iter);
        foundEqual = true;
        break;
      }
    }

    if(!foundEqual) {
      groups.emplace_back(
        std::vector<T> {*iter}
      );
    }
  }

  return groups;
}

} // namespace temple
} // namespace Scine

#endif
