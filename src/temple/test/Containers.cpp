/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "temple/Stringify.h"
#include "temple/constexpr/JSF.h"
#include "temple/OrderedPair.h"
#include "temple/Poset.h"
#include "temple/Random.h"

#include <algorithm>

extern temple::jsf::Generator<> generator;

BOOST_AUTO_TEST_CASE(OrderedPairTests) {
  temple::OrderedPair<unsigned> a {14u, 3u};
  BOOST_CHECK(a.front() < a.back());
  BOOST_CHECK(a.first < a.second);
  BOOST_CHECK(
    std::is_sorted(std::begin(a), std::end(a))
  );
}

BOOST_AUTO_TEST_CASE(PosetTests) {
  constexpr unsigned nTests = 10;
  constexpr unsigned nValues = 20;
  constexpr unsigned maxValue = 40;

  auto compareFirstDigit = [](const unsigned a, const unsigned b) -> bool {
    return (a / 10) < (b / 10);
  };

  std::vector<unsigned> unorderedValues;
  for(unsigned testNum = 0; testNum < nTests; ++testNum) {
    unorderedValues = temple::random::getN<unsigned>(0, maxValue, nValues, generator.engine);
    temple::Poset<unsigned> f {unorderedValues};

    f.orderUnordered(compareFirstDigit);
    f.orderUnordered(std::less<>());
    f.finalize();

    BOOST_CHECK(
      std::is_sorted(
        std::begin(f),
        std::end(f),
        [](const auto& a, const auto& b) -> bool {
          return std::less<>()(
            a.values.front(),
            b.values.front()
          );
        }
      )
    );
  }
}
