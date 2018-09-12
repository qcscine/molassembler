#include <boost/test/unit_test.hpp>

#include "temple/Cache.h"

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
  mutable temple::Cache<std::string> _cache {
    std::make_pair(
      "bigNumber",
      [this]() {
        return (this -> _determineMe());
      }
    )
  };

  unsigned _determineMe() {
    return ackermann(3, 2);
  }

public:
  unsigned getAckermann() const {
    /* 4 */
    return _cache.getGeneratable<unsigned>("bigNumber");
  }

  void changeCacheValue() {
    /* 5 */
    _cache.changeGeneratable<unsigned>(
      "bigNumber",
      [](unsigned* value) {
        *value = 4;
      }
    );
  }
};

BOOST_AUTO_TEST_CASE( cache_all ) {
  using namespace std::string_literals;

  /* 1 */
  temple::Cache<std::string> cache;

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
