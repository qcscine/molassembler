#include <boost/test/unit_test.hpp>

#include "temple/Adaptors.h"
#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"

#define INCLUDE_INVOKE_WITH_BOOST_TUPLE
#include "temple/Invoke.h"

#include <iostream>

template<class Container, class UnaryPredicate>
bool trial_all_of(const Container& container, UnaryPredicate&& predicate) PURITY_WEAK;

template<class Container, class UnaryPredicate>
bool trial_all_of(const Container& container, UnaryPredicate&& predicate) {
  for(const auto& element : container) {
    if(!temple::invoke(predicate, element)) {
      return false;
    }
  }

  return true;
}

unsigned weirdIntDiv(unsigned a, unsigned b) PURITY_STRONG;
unsigned weirdIntDiv(unsigned a, unsigned b) {
  return a / (b + 1);
}

unsigned plusFive(unsigned a) PURITY_STRONG;
unsigned plusFive(unsigned a) {
  return a + 5;
}

BOOST_AUTO_TEST_CASE(adaptors) {
  using namespace temple;

  std::vector<unsigned> a {4, 1, 3, 5};
  std::vector<unsigned> b {5, 6, 7, 8};

  auto loudCompare = [](unsigned a, unsigned b) -> bool {
    std::cout << a << " < " << b << ": " << (a < b) << "\n";
    return (a < b);
  };

  bool pairwiseSmaller = all_of(
    zip(a, b),
    make_tuple_callable(loudCompare)
  );

  bool bSortedAsc = all_of(
    sequentialPairs(b),
    make_tuple_callable(std::less<>())
  );

  BOOST_CHECK(pairwiseSmaller && bSortedAsc);

  BOOST_CHECK(
    trial_all_of(
      zip(a, b),
      std::less<>()
    )
  );

  BOOST_CHECK(
    sum(
      transform(
        allPairs(a),
        std::plus<>()
      )
    ) == 5u + 7u + 9u + 4u + 6u + 8u
  );

  BOOST_CHECK(
    sum(
      transform(
        allPairs(a),
        weirdIntDiv
      )
    ) == 3u
  );

  BOOST_CHECK(
    sum(
      transform(a, plusFive)
    ) == 33u
  );

  BOOST_CHECK(
    trial_all_of(
      sequentialPairs(b),
      std::less<>()
    )
  );

  BOOST_CHECK(
    trial_all_of(
      a,
      [](unsigned a) {
        return a < 6;
      }
    )
  );

  BOOST_CHECK(
    trial_all_of(
      a,
      std::bind(std::less<>(), std::placeholders::_1, 6u)
    )
  );
}
