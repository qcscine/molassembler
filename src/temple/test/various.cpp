#include <boost/test/results_collector.hpp>
#include <boost/test/unit_test.hpp>

#include "boost/optional.hpp"

#include "temple/Functional.h"
#include "temple/OperatorSuppliers.h"
#include "temple/Stringify.h"
#include "temple/constexpr/Numeric.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

template<typename T>
std::ostream& operator << (std::ostream& os, const boost::optional<T>& valueOptional) {
  if(valueOptional) {
    os << "Some " << valueOptional.value();
  } else {
    os << "None";
  }

  return os;
}

template<typename NullaryCallable>
double timeNullaryCallable(
  NullaryCallable&& function
) {
  using namespace std::chrono;

  time_point<system_clock> start, end;
  std::array<unsigned, 1000> timings;
  for(unsigned N = 0; N < 1000; ++N) {
    start = system_clock::now();

    function();

    end = system_clock::now();
    duration<double> elapsed = end - start;

    timings.at(N) = elapsed.count() * 1e9;
  }

  return temple::average(timings);
}


BOOST_AUTO_TEST_CASE(setRemovalDefect) {
  /* Sets and maps cannot use std::remove_if! Not a defect. */

  auto testSet = temple::copy_if(
    std::set<unsigned> {5, 2, 9, 1, 3, 4, 8},
    [](const auto& value) -> bool {
      return value % 2 != 0;
    }
  );

  BOOST_CHECK((testSet == std::set<unsigned> {1, 3, 5, 9}));
}

BOOST_AUTO_TEST_CASE(selectTestCases) {
  auto ragged2D = std::vector<
    std::vector<unsigned>
  > {
    {4, 9, 3},
    {2, 11, 6},
    {30, 4},
    {9, 44, 33, 12}
  };

  auto smallestVector = temple::select(
    ragged2D,
    std::less<>(),
    [](const auto& group) -> unsigned {
      return group.size();
    }
  );

  BOOST_CHECK((*smallestVector == std::vector<unsigned> {30, 4}));

  auto largestVector = temple::select(
    ragged2D,
    std::greater<>(),
    [](const auto& group) -> unsigned {
      return group.size();
    }
  );

  BOOST_CHECK((*largestVector == std::vector<unsigned> {9, 44, 33, 12}));
}

BOOST_AUTO_TEST_CASE(stringifyTests) {
  std::vector<
    std::map<
      unsigned,
      std::pair<
        int,
        double
      >
    >
  > complicatedStructure {
    {
      {
        4,
        {-2, 4.5}
      },
      {
        9,
        {-10, 0.1}
      }
    },
    {
      {
        100,
        {5, 1.2}
      }
    }
  };

  std::cout << temple::stringify(complicatedStructure) << std::endl;

  enum class SampleEnum {
    SomeValue,
    AnotherValue
  };

  std::tuple<
    std::vector<unsigned>,
    int,
    std::pair<double, int>,
    boost::optional<unsigned>,
    SampleEnum
  > someTuple {
    {0, 4, 13, 4},
    -4,
    {0.9, -100},
    boost::optional<unsigned> {10},
    SampleEnum::SomeValue
  };

  std::cout << temple::stringify(someTuple) << std::endl;
}

struct TupleLike : temple::crtp::AllOperatorsFromTupleMethod<TupleLike> {
  unsigned x;
  int f;

  TupleLike(unsigned x_, int f_) : x(x_), f(f_) {}

  auto tuple() const {
    return std::tie(x, f);
  }
};

BOOST_AUTO_TEST_CASE(crtpTests) {
  TupleLike a {4u, -3}, b {9u, 4}, c {9u, 4};
  BOOST_CHECK(b == c);
  BOOST_CHECK(a < b);
  BOOST_CHECK(b > a);
  BOOST_CHECK(b <= c);
  BOOST_CHECK(b >= c);
}
