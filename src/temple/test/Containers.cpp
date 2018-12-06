/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "temple/OrderedPair.h"

#include <algorithm>

BOOST_AUTO_TEST_CASE(OrderedPairTests) {
  temple::OrderedPair<unsigned> a {14u, 3u};
  BOOST_CHECK(a.front() < a.back());
  BOOST_CHECK(a.first < a.second);
  BOOST_CHECK(
    std::is_sorted(std::begin(a), std::end(a))
  );
}
