/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Abstractions of variable-depth loops
 */

#ifndef INCLUDE_TEMPLE_LOOPS_H
#define INCLUDE_TEMPLE_LOOPS_H

#include <vector>
#include <algorithm>
#include <stdexcept>

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Loops {
namespace Detail {

template<typename T, typename F, typename Predicate>
void descend(
  F&& f,
  Predicate&& p,
  const unsigned depth,
  const T start,
  const T end,
  std::vector<T>& indices
) {
  if(depth == 0) {
    f(indices);
    return;
  }

  const unsigned I = indices.size() - depth;
  for(indices[I] = start; indices[I] < end; ++indices[I]) {
    if(p(std::begin(indices), std::begin(indices) + I, indices[I])) {
      descend(f, p, depth - 1, start, end, indices);
    }
  }
}

} // namespace Detail

/**
 * @brief Call a function for every distinct combination of N indices within a range
 *
 * An abstraction of i != j != k != ... with {i, j, k, ...} in [start, end),
 * calls @p f with a vector of size @p count distinct indices.
 */
template<typename T, typename F>
void different(F&& f, const unsigned count, const T start, const T end) {
  auto descentPredicate = [](auto first, auto last, const T i) -> bool {
    return std::find(first, last, i) == last;
  };

  if(count == 0) {
    throw std::logic_error("Zero loop indices makes no sense");
  }

  if(start > end) {
    throw std::logic_error("Start is bigger than end, you asked for an infinite loop");
  }

  std::vector<T> indices(count, start);
  for(; indices[0] < end; ++indices[0]) {
    Detail::descend(f, descentPredicate, count - 1, start, end, indices);
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

} // namespace Loops
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
