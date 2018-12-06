/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "temple/Variadic.h"

#include <set>
#include <vector>

BOOST_AUTO_TEST_CASE(concatenateTests) {
  std::set<unsigned> f {5, 9, 3}; // Note these get reordered
  std::vector<unsigned> h {9, 7, 4};

  BOOST_CHECK(temple::variadic::sizeLowerBound(f, h) == 6u);
  auto concatenated = temple::variadic::concatenate(f, h);
  BOOST_CHECK((concatenated == std::vector<unsigned> {3, 5, 9, 9, 7, 4}));
}
