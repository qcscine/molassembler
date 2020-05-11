/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Partitioner.h"

#include "molassembler/Temple/Functional.h"

#include <cassert>

namespace Scine {
namespace Molassembler {
namespace Shapes {

Partitioner::Partitioner(const unsigned s, const unsigned e) : S(s), E(e), mapping(S * E) {
  assert(s != 0 && e != 0);
  const unsigned numElements = S * E;
  for(unsigned i = 0; i < numElements; ++i) {
    mapping[i] = i / E;
  }

  assert(isOrderedMapping(mapping));
}

bool Partitioner::next_partition() {
  if(mapping.size() == 1) {
    return false;
  }

  assert(mapping.size() > 1);
  assert(mapping.front() == 0);
  assert(mapping.size() == S * E);
  std::vector<unsigned> counts(S, E);

  /* The last position is not free, it is predetermined by all other counts.
   * Decrement the count of what is set there and move left.
   */
  auto it = std::end(mapping) - 1;
  assert(*it < S);
  --counts[*it];
  --it;

  /* Find an incrementable position from the back up to and including the
   * second position. The first position is not free to choose, it must be 0.
   */
  bool incremented = false;
  for(const auto first = std::begin(mapping); it != first; --it) {
    // Decrement the count at this position
    assert(*it < S);
    --counts[*it];

    /* An incrementable position can choose a successor larger than its current
     * number. A successor is only possible if all preceding numbers have been
     * placed at least once.
     */
    bool predecessorsPlaced = true;
    for(unsigned i = 0; i < *it; ++i) {
      if(counts[i] == 0) {
        predecessorsPlaced = false;
        break;
      }
    }

    if(predecessorsPlaced) {
      // We know counts(0, ..., (*it - 1)) > 0
      for(unsigned i = *it + 1; i < S; ++i) {
        /* A successor i can be chosen if
         * - all i-1, i-2, ..., 2 have been chosen before, i.e. all of the
         *   preceding counts to i are > 0. If any of them is 0, then neither i
         *   nor its increment could be chosen, so break.
         * - it has been chosen less than E times, i.e. its count is < E
         */
        assert(i != 0);
        if(counts[i - 1] == 0) {
          break;
        }

        if(counts[i] < E) {
          // We can increment here! Do so:
          incremented = true;
          *it = i;
          ++counts[i];
          break;
        }
      }

      if(incremented) {
        /* This is not a loop condition because then the iterator it would have
         * been needlessly decremented further.
         */
        break;
      }
    }

    // No incrementable successor found, keep looking at the previous position
  }

  if(!incremented) {
    return false;
  }

  /* Choose the lowest partition forward from it */
  ++it;
  for(const auto end = std::end(mapping); it != end; ++it) {
    for(unsigned i = 0; i < S; ++i) {
      if(counts[i] < E) {
        *it = i;
        ++counts[i];
        break;
      }
    }
  }

  return true;
}

std::vector<
  std::vector<unsigned>
> Partitioner::partitions() const {
  std::vector<
    std::vector<unsigned>
  > groups(S);
  for(auto& group : groups) {
    group.reserve(E);
  }

  const unsigned numElements = S * E;
  for(unsigned i = 0; i < numElements; ++i) {
    assert(mapping[i] < S);
    groups[mapping[i]].push_back(i);
  }

  assert(Temple::all_of(groups, [&](const auto& g) { return g.size() == E; }));

  return groups;
}

bool Partitioner::isOrderedMapping(const std::vector<unsigned>& mapping) {
  /* We only need to check lexicographic ordering of first position. The
   * iterators CANNOT compare equal (since they are looking for different
   * indices or mapping has been setup incorrectly) so unless any of them are
   * pairwise bigger all of them are sequentially smaller indicating
   * lexicographic ordering.
   *
   * The remaining positions' order does not matter as the first position of a
   * group index fully indicates lexicographic ordering.
   */
  // The first position is never altered, so we can skip it in checking
  assert(mapping[0] == 0);

  const unsigned N = mapping.size();
  unsigned boundary = 1;
  for(unsigned i = 1; i < N; ++i) {
    const unsigned x = mapping[i];
    if(x > boundary) {
      return false;
    }

    if(x == boundary) {
      ++boundary;
    }
  }

  return true;
}

} // namespace Shapes
} // namespace Molassembler
} // namespace Scine
