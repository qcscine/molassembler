#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_SET_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_SET_ALGORITHMS_H

#include <algorithm>

/*! @file
 *
 * @brief Shortcuts to call STL set algorithms.
 */

namespace temple {

template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> set_intersection(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
) {
  std::set<T, Comparator<T>, Allocator<T>> returnSet;

  std::set_intersection(
    std::begin(a),
    std::end(a),
    std::begin(b),
    std::end(b),
    std::inserter(returnSet, std::end(returnSet)),
    Comparator<T>()
  );

  return returnSet;
}

template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> set_union(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
) {
  std::set<T, Comparator<T>, Allocator<T>> returnSet;

  std::set_union(
    std::begin(a),
    std::end(a),
    std::begin(b),
    std::end(b),
    std::inserter(returnSet, std::end(returnSet)),
    Comparator<T>()
  );

  return returnSet;
}

template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> set_symmetric_difference(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
) {
  std::set<T, Comparator<T>, Allocator<T>> returnSet;

  std::set_symmetric_difference(
    std::begin(a),
    std::end(a),
    std::begin(b),
    std::end(b),
    std::inserter(returnSet, std::end(returnSet)),
    Comparator<T>()
  );

  return returnSet;
}

} // namespace temple

#endif
