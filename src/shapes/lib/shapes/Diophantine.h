/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Tools to enumerate all solutions of simple diophantine equations
 */

#include <vector>

namespace Scine {
namespace shapes {

/*! @brief Functions to help treat a particular linear diophantine equation
 *
 * Here we want to treat a special case of the linear diophantine equation
 * @math{\sum_i a_i x_i = b} with @math{a_i > 0} and @math{x_i >= 0}.
 */
namespace diophantine {

/*! @brief Checks whether a diophantine has a solution
 *
 * @param a The constant coefficients, in no particular order and without
 *   value constraints. The list may not be empty, though.
 * @param b The sought inner product result.
 *
 * @returns whether `gcd(a_1, ..., a_n)` is a divisor of b. If so, the diophantine
 *   has a solution.
 */
bool has_solution(
  const std::vector<unsigned>& a,
  const int b
);

/*! @brief Finds the next solution of the diophantine equation, if it exists
 *
 * Finds the next solution of the linear diophantine @math{\sum_i a_ix_i = b} with
 * positive @math{a_i}, non-negative @math{x_i} and positive @math{b}.
 *
 * Use the following pattern:
 * @code{cpp}
 * std::vector<unsigned> x;
 * const std::vector<unsigned> a {4, 3, 2};
 * const int b = 12;
 *
 * if(first_solution(x, a, b)) {
 *   do {
 *     // Do something with x
 *   } while(next_solution(x, a, b));
 * }
 * @endcode
 *
 * @warning Solutions to diophantines, even particularly easy ones like this
 * constrained, linear diophantine, are generally complex. The approach taken
 * here is more or less brute-force enumeration and does not scale well for
 * many coefficients `a`! See https://arxiv.org/pdf/math/0010134.pdf for a
 * sketch on how this could be properly solved.
 *
 * @param x Filled coefficient vector of equal size as a, whose values
 *   represent a solution to the diophantine (i.e. `inner_product(a, x) == b`)
 * @param a The constant coefficients, strictly descending (no duplicate
 *   values) and non-zero.
 * @param b The sought inner product result
 */
bool next_solution(
  std::vector<unsigned>& x,
  const std::vector<unsigned>& a,
  const int b
);

/*! @brief Finds the first solution of the diophantine equation, if it exists
 *
 * Finds the first solution of the linear diophantine @math{\sum_i a_ix_i = b} with
 * positive @math{a_i}, non-negative @math{x_i} and positive @math{b}.
 *
 * Use the following pattern:
 * @code{cpp}
 * std::vector<unsigned> x;
 * const std::vector<unsigned> a {4, 3, 2};
 * const int b = 12;
 *
 * if(first_solution(x, a, b)) {
 *   do {
 *     // Do something with x
 *   } while(next_solution(x, a, b));
 * }
 * @endcode
 *
 * @warning Solutions to diophantines, even particularly easy ones like this
 * constrained, linear diophantine, are generally complex. The approach taken
 * here is more or less brute-force enumeration and does not scale well for
 * many coefficients @p `a`! See https://arxiv.org/pdf/math/0010134.pdf for a
 * sketch on how this could be properly solved.
 *
 * @param x Vector to store the first solution into (resize and fill is handled
 *   by this function)
 * @param a The constant coefficients, strictly descending (no duplicate
 *   values) and non-zero.
 * @param b The sought inner product result
 */
bool first_solution(
  std::vector<unsigned>& x,
  const std::vector<unsigned>& a,
  const int b
);

} // namespace diophantine
} // namespace shapes
} // namespace Scine
