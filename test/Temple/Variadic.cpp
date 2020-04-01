/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Temple/Variadic.h"

#include <set>
#include <vector>

using namespace Scine;

BOOST_AUTO_TEST_CASE(concatenateTests) {
  std::set<unsigned> f {5, 9, 3}; // Note these get reordered
  std::vector<unsigned> h {9, 7, 4};

  BOOST_CHECK(temple::variadic::sizeLowerBound(f, h) == 6u);
  auto concatenated = temple::variadic::concatenate(f, h);
  BOOST_CHECK((concatenated == std::vector<unsigned> {3, 5, 9, 9, 7, 4}));
}
