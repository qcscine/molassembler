/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details. for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Shapes/Diophantine.h"

using namespace Scine::shapes;

BOOST_AUTO_TEST_CASE(Diophantine) {
  std::vector<unsigned> x;
  const std::vector<unsigned> a {4, 3, 2};
  const int b = 12;

  const std::vector<
    std::vector<unsigned>
  > expectedX {
    {0, 0, 6},
    {0, 2, 3},
    {0, 4, 0},
    {1, 0, 4},
    {1, 2, 1},
    {2, 0, 2},
    {3, 0, 0}
  };

  BOOST_REQUIRE(diophantine::first_solution(x, a, b));
  unsigned i = 0;
  do {
    BOOST_CHECK(x == expectedX.at(i));
    ++i;
  } while(diophantine::next_solution(x, a, b));
  BOOST_REQUIRE_EQUAL(i, expectedX.size());
  BOOST_REQUIRE(x == std::vector<unsigned> (3, 0));
}
