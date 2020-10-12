/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Variadic.h"

#include <set>
#include <vector>

using namespace Scine::Molassembler;

BOOST_AUTO_TEST_CASE(ConcatenateTests, *boost::unit_test::label("Temple")) {
  std::set<unsigned> f {5, 9, 3}; // Note these get reordered
  std::vector<unsigned> h {9, 7, 4};

  BOOST_CHECK(Temple::variadic::sizeLowerBound(f, h) == 6u);
  auto concatenated = Temple::variadic::concatenate(f, h);
  BOOST_CHECK((concatenated == std::vector<unsigned> {3, 5, 9, 9, 7, 4}));
}
