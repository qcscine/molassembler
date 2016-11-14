#include <iostream>
#include <vector>
#include <map>
#include <typeinfo>
#include <chrono>

#include "Cache.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

// Some help functions for testing
template<typename T>
unsigned time_wrap_nullary_callable(
    const T& nullary_callable
) {
    using namespace std::chrono;

    time_point<system_clock> start, end;
    start = system_clock::now();

    nullary_callable();
    
    end = system_clock::now();
    duration<double> elapsed = end - start;

    return elapsed.count() * 1e3;
}

unsigned int ackermann(unsigned int m, unsigned int n) {
  if (m == 0) {
    return n + 1;
  }
  if (n == 0) {
    return ackermann(m - 1, 1);
  }
  return ackermann(m - 1, ackermann(m, n - 1));
}


/* Cache members to test
 *
 * t # name
 * y 1 default constructor
 * y 2 init list constructor
 * y 3 add
 * y 4 getGeneratable
 * y 5 changeGeneratable
 * y 6 invalidate (blanket)
 * y 7 invalidate (selective)
 * y 8 getOption
 * y 9 has
 *
 */

/* Sample class with mutable Cache member and a generatable cache object
 * Also includes an example of how to modify a cache object
 */
class Foo {
private:
  /* 2 */
  mutable Cache _cache {
    std::make_pair(
      "bigNumber",
      [this]() {
        return (this -> _determineMe());
      }
    )
  };

  unsigned _determineMe() {
    return ackermann(4, 1);
  }

public:
  unsigned getAckermann41() const {
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

using namespace std::string_literals;

BOOST_AUTO_TEST_CASE( cache_all ) {
  /* 1 */
  Cache cache;

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
  std::cout << "Calculation: " << time_wrap_nullary_callable(
    [&]() {
      bar.getAckermann41();
    }
  ) << "ms" << std::endl;
  std::cout << "From cache: " << time_wrap_nullary_callable(
    [&]() {
      bar.getAckermann41();
    }
  ) << "ms" << std::endl;

  bar.changeCacheValue();
  BOOST_CHECK(bar.getAckermann41() == 4);
}
