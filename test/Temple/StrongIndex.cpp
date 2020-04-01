/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Temple/StrongIndex.h"
#include <vector>

using namespace Scine;

BOOST_AUTO_TEST_CASE(StrongIndices) {
  struct foo_tag;
  using Foo = temple::StrongIndex<foo_tag, unsigned>;

  struct bar_tag;
  using Bar = temple::StrongIndex<bar_tag, unsigned>;

  Foo f(4);
  ++f;
  BOOST_CHECK(f == 5u);

  Foo g(5);
  BOOST_CHECK(f == g);

  // This is an eyesore, possible because of non-explicit conversion
  Bar h(5);
  BOOST_CHECK(g == h);

  // This should also not be permitted
  g = h;
}
