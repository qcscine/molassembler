/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Shapes/Diophantine.h"

using namespace Scine::Molassembler::Shapes;

BOOST_AUTO_TEST_CASE(DiophantineExample, *boost::unit_test::label("Shapes")) {
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

  BOOST_REQUIRE(Diophantine::first_solution(x, a, b));
  unsigned i = 0;
  do {
    BOOST_CHECK(x == expectedX.at(i));
    ++i;
  } while(Diophantine::next_solution(x, a, b));
  BOOST_REQUIRE_EQUAL(i, expectedX.size());
  BOOST_REQUIRE(x == std::vector<unsigned> (3, 0));
}
