/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Abstractions of variable-depth loops
 */

#ifndef INCLUDE_TEMPLE_LOOPS_H
#define INCLUDE_TEMPLE_LOOPS_H

#include <vector>
#include <algorithm>
#include <stdexcept>

namespace Scine {
namespace temple {
namespace loops {
namespace detail {

template<typename T, typename F>
void different(F&& f, const unsigned depth, const T start, const T end, std::vector<T>& indices) {
  if(depth == 0) {
    f(indices);
    return;
  } else {
    const unsigned I = indices.size() - depth;
    for(indices[I] = start; indices[I] < end; ++indices[I]) {
      if(std::find(std::begin(indices), std::begin(indices) + I, indices[I]) != std::begin(indices) + I) {
        continue;
      }

      different(f, depth - 1, start, end, indices);
    }
  }
}

} // namespace detail

/**
 * @brief Call a function for every distinct combination of N indices within a range
 *
 * An abstraction of i != j != k != ... with {i, j, k, ...} in [start, end),
 * calls @p f with a vector of size @p count distinct indices.
 */
template<typename T, typename F>
void different(F&& f, const unsigned count, const T start, const T end) {
  if(count == 0) {
    throw std::logic_error("Zero loop indices makes no sense");
  }

  if(start > end) {
    throw std::logic_error("Start is bigger than end, you asked for an infinite loop");
  }

  std::vector<T> indices(count, start);
  for(; indices[0] < end; ++indices[0]) {
    detail::different(f, count - 1, start, end, indices);
  }
}

/**
 * @brief Call a function for every distinct combination of N indices within a range
 *
 * An abstraction of i != j != k != ... with {i, j, k, ...} in [0, end),
 * calls @p f with a vector of size @p count distinct indices.
 */
template<typename T, typename F>
void different(F&& f, const unsigned count, const T end) {
  different(std::forward<F>(f), count, T {0}, end);
}

} // namespace loops
} // namespace temple
} // namespace Scine

#endif
