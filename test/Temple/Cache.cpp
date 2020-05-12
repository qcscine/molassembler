/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Cache.h"

using namespace Scine::Molassembler;

unsigned int ackermann(unsigned int m, unsigned int n) {
  if (m == 0) {
    return n + 1;
  }
  if (n == 0) {
    return ackermann(m - 1, 1);
  }
  return ackermann(m - 1, ackermann(m, n - 1));
}

/* Sample class with mutable Cache member and a generatable cache object
 * Also includes an example of how to modify a cache object
 */
class Foo {
private:
  /* 2 */
  mutable Temple::Cache<std::string> cache_ {
    std::make_pair(
      "bigNumber",
      [this]() {
        return (this -> determineMe_());
      }
    )
  };

  unsigned determineMe_() {
    return ackermann(3, 2);
  }

public:
  unsigned getAckermann() const {
    /* 4 */
    return cache_.getGeneratable<unsigned>("bigNumber");
  }

  void changeCacheValue() {
    /* 5 */
    cache_.changeGeneratable<unsigned>(
      "bigNumber",
      [](unsigned* value) {
        *value = 4;
      }
    );
  }
};

BOOST_AUTO_TEST_CASE(SimpleCacheTest, *boost::unit_test::label("Temple")) {
  using namespace std::string_literals;

  /* 1 */
  Temple::Cache<std::string> cache;

  /* 3 */
  std::vector<std::string> keys {"number", "string", "vector"};
  cache.add("number", 5);
  cache.add("string", "fsldkf"s);
  cache.add("vector", std::vector<unsigned>({4, 9}));

  /* 8 */
  if(auto number = cache.getOption<int>("number")) { // op (bool) is true
    BOOST_CHECK(number.value() == 5);
  }

  if(auto number = cache.getOption<int>("non-existent number")) { // op (bool) is false
    // this should not be executed
    BOOST_REQUIRE(false);
  }

  /* 9 */
  BOOST_CHECK(
    std::all_of(
      keys.begin(),
      keys.end(),
      [&cache](const std::string& key) {
        return cache.has(key);
      }
    )
  );

  /* 6 */
  cache.invalidate("number");
  BOOST_CHECK(!cache.has("number"));

  /* 7 */
  cache.invalidate();
  BOOST_CHECK(
    std::none_of(
      keys.begin(),
      keys.end(),
      [&cache](const std::string& key) {
        return cache.has(key);
      }
    )
  );

  /* 2, 4, 5 */
  Foo bar;
  bar.getAckermann();

  // test modification of the cache
  bar.changeCacheValue();
  BOOST_CHECK(bar.getAckermann() == 4);
}
