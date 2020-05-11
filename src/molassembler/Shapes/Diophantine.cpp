/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Diophantine.h"

#include "boost/integer/common_factor.hpp"
#include <functional>
#include <algorithm>
#include <numeric>
#include <cassert>

namespace Scine {
namespace Shapes {
namespace Diophantine {

bool has_solution(
  const std::vector<unsigned>& a,
  const int b
) {
  // Calculate gcd for all elements in a
  assert(!a.empty());
  unsigned collectiveGCD = a.front();
  for(auto it = ++std::begin(a); it != std::end(a); ++it) {
    collectiveGCD = boost::integer::gcd(collectiveGCD, *it);
  }

  /* For the diophantine to have a solution, the collective GCD has to evenly
   * divide b
   */
  return b % collectiveGCD == 0;
}

bool next_solution(
  std::vector<unsigned>& x,
  const std::vector<unsigned>& a,
  const int b
) {
  // a must be ordered descending
  assert(
    std::is_sorted(
      std::begin(a),
      std::end(a),
      std::greater<>()
    )
  );
  // a and x sizes must match
  assert(a.size() == x.size());
  assert(!x.empty());
  assert(b > 0);

  if(x.size() == 1) {
    return false;
  }

  assert(x.size() >= 2);
  auto q = [&]() -> int {
    return std::inner_product(
      std::begin(a),
      --std::end(a),
      std::begin(x),
      0
    );
  };

  /* To find the next solution, we increment the left n - 1 columns of a
   * through all combinations giving inner products <= b until the difference
   * is divisible by the last a.
   */
  const auto reverseStart = ++std::rbegin(x);
  const auto reverseEnd = std::rend(x);
  assert(reverseStart != reverseEnd);
  int Q = 0;
  do {
    // Increment the left n - 1 columns of x once!
    auto it = reverseStart;
    assert(it != reverseEnd);
    for(; it != reverseEnd; ++it) {
      ++(*it);

      Q = q();
      if(Q <= b) {
        break;
      }

      *it = 0;
    }

    // Check to make sure we are not at the last solution
    if(it == reverseEnd) {
      // Past last solution!
      x.back() = 0;
      return false;
    }

  } while(Q > b || (b - Q) % a.back() != 0);

  x.back() = (b - Q) / a.back();
  return true;
}

bool first_solution(
  std::vector<unsigned>& x,
  const std::vector<unsigned>& a,
  const int b
) {
  assert(!a.empty());

  if(a.size() == 1) {
    x = {b / a.front()};
    return true;
  }

  assert(a.size() >= 2);
  x.resize(a.size());
  std::fill(
    std::begin(x),
    std::end(x),
    0
  );

  if(b % a.back() == 0) {
    x.back() = b / a.back();
    return true;
  }

  return next_solution(x, a, b);
}

} // namespace Diophantine
} // namespace Shapes
} // namespace Scine
