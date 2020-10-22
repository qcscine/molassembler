/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/StrongIndex.h"
#include <vector>

using namespace Scine::Molassembler;

BOOST_AUTO_TEST_CASE(StrongIndices, *boost::unit_test::label("Temple")) {
  struct FooTag;
  using Foo = Temple::StrongIndex<FooTag, unsigned>;

  struct BarTag;
  using Bar = Temple::StrongIndex<BarTag, unsigned>;

  Foo f(4);
  ++f;
  BOOST_CHECK(f == 5);

  Foo g(5);
  BOOST_CHECK(f == g);

  // This is an eyesore, possible because of non-explicit conversion
  Bar h(5);
  BOOST_CHECK(g == h);

  // This should also not be permitted
  g = h;
}
