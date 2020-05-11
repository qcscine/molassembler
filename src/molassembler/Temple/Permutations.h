/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides functionality related to permutations.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_PERMUTATIONS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_PERMUTATIONS_H

#include <algorithm>
#include <functional>

namespace Scine {
namespace Molassembler {
namespace Temple {

/*! @brief Calculate the index of permutation of elements in a container
 *
 * @complexity{@math{\Theta(N^2)}}
 * @note Requires that Container implements operator[](U), where U is implicitly
 *   convertible from std::size_t.
 */
template<class Container>
std::size_t permutationIndex(const Container& container) {
  const std::size_t size = container.size();

  std::size_t index = 0;
  std::size_t position = 2;// position 1 is paired with factor 0 and so is skipped
  std::size_t factor = 1;

  for(std::size_t p = size - 2; p != std::numeric_limits<std::size_t>::max(); --p) {
    std::size_t largerSuccessors = 0;

    for(std::size_t q = p + 1; q < size; ++q) {
      if(container[q] < container[p]) {
        ++largerSuccessors;
      }
    }

    index += (largerSuccessors * factor);
    factor *= position;
    ++position;
  }

  return index;
}

namespace InPlace {

//! Calls std::next_permutation
template<class Container>
bool next_permutation(Container& container) {
  return std::next_permutation(
    std::begin(container),
    std::end(container)
  );
}

//! Calls std::prev_permutation
template<class Container>
bool prev_permutation(Container& container) {
  return std::prev_permutation(
    std::begin(container),
    std::end(container)
  );
}

/*!
 * @brief For when you have to implement variable-depth for loops, each with
 *   different limits.
 *
 * Increments a container containing indices for each depth according to a const
 * container containing the limits for each. Limits are exclusive!
 *
 * E.g. for the limits 232, increments through the sequence
 * 000 -> 001 -> 010 -> 011 -> 020 -> 021 -> 100 -> 101 -> 110 -> 111 -> 120 ...
 *
 * @note Requires that Container implements operator[](U), where U is implicitly
 *   convertible from unsigned.
 *
 * @complexity{@math{\Theta(N)}}
 */
template<class Container>
bool nextCombinationPermutation(
  Container& toPermute,
  const Container& limits
) {
  assert(toPermute.size() == limits.size());
  const unsigned cols = toPermute.size();

  // Check if all columns are full
  bool allFull = true;
  for(unsigned i = 0; i < cols; ++i) {
    if(toPermute[i] != limits[i]) {
      allFull = false;
      break;
    }
  }

  if(allFull) {
    return false;
  }

  // Make next permutation
  for(int i = cols - 1; i >= 0; --i) {
    if(toPermute[i] == limits[i]) {
      toPermute[i] = 0;
    } else {
      ++toPermute[i];
      return true;
    }
  }

  return true;
}

} // namespace InPlace
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
