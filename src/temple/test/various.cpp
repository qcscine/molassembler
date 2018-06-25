#define BOOST_TEST_MODULE VariousTempleTestsModule
#define BOOST_TEST_DYN_LINK

#include <boost/test/results_collector.hpp>
#include <boost/test/unit_test.hpp>

#include "boost/optional.hpp"

#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Enumerate.h"
#include "temple/MemberFetcher.h"
#include "temple/Random.h"
#include "temple/VectorView.h"
#include "temple/Stringify.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

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


BOOST_AUTO_TEST_CASE( enumerateTests) {
  std::vector<unsigned> testVec {5, 2, 3, 4};

  bool pass = true;
  for(const auto& enumPair : enumerate(testVec)) {
    if(testVec.at(enumPair.index) != enumPair.value) {
      pass = false;
      break;
    }
  }

  BOOST_CHECK(pass);

  auto weirdSum = temple::sum(
    temple::mapToVector(
      enumerate(testVec),
      [](const auto& enumPair) -> unsigned {
        return enumPair.index + enumPair.value;
      }
    )
  );

  BOOST_CHECK(weirdSum == 5 + 3 + 5 + 7);
}

BOOST_AUTO_TEST_CASE(memberFetcherTests) {
  std::vector<
    std::vector<unsigned>
  > testClass {
    {4, 3, 9, 10},
    {1, 2, 3, 4, 5},
    {11},
    {44, 93}
  };

  auto memberFetcherClass = temple::getMember(
    testClass,
    [](const auto& container) -> unsigned {
      return container.size();
    }
  );

  BOOST_CHECK_MESSAGE(
    temple::sum(memberFetcherClass) == 12u,
    "Member fetcher does not work as expected, begin()-end() contains {"
      << temple::condenseIterable(memberFetcherClass)
      << "} instead of expected {4, 5, 1, 2}"
  );

  BOOST_CHECK_MESSAGE(
    temple::min(memberFetcherClass) == 1
    && temple::max(memberFetcherClass) == 5,
    "Usage of min/max does not return expected result"
  );

  // Is this approach more efficient than copying member information?
  double copyingTime = timeNullaryCallable(
    [&]() {
      auto mappedData = temple::map(
        testClass,
        [](const auto& subVector) -> unsigned {
          return subVector.size();
        }
      );

      for(const auto& size : mappedData) {
        static_cast<void>(size);
      }
    }
  );

  double memberFetcherTime = timeNullaryCallable(
    [&]() {
      auto fetcher = temple::getMember(
        testClass,
        [](const auto& subVector) -> unsigned {
          return subVector.size();
        }
      );

      for(const auto& size : fetcher) {
        static_cast<void>(size);
      }
    }
  );

  BOOST_CHECK_MESSAGE(
    memberFetcherTime < copyingTime,
    "Copying is faster than fetching members: copy = " << copyingTime << " vs. "
    << memberFetcherTime << " member fetch "
  );

  // Interoperability with VectorView

  auto filteredView = temple::view_filter(
    testClass,
    [](const auto& subVector) -> bool {
      return subVector.size() < 3;
    }
  );

  auto minSize = temple::min(
    temple::getMember(
      filteredView,
      [](const auto& subVector) -> unsigned {
        return subVector.size();
      }
    )
  );

  BOOST_CHECK_MESSAGE(
    minSize == 4,
    "Min size of filteredView returns " << minSize << ", not 4"
  );

  auto maxSize = temple::max(
    temple::getMember(
      filteredView,
      [](const auto& subVector) -> unsigned {
        return subVector.size();
      }
    )
  );

  BOOST_CHECK_MESSAGE(
    maxSize == 5,
    "Max size of filteredView returns " << maxSize << ", not 5"
  );
}

BOOST_AUTO_TEST_CASE(setRemovalDefect) {
  /* Sets and maps cannot use std::remove_if! Not a defect. */

  auto testSet = temple::moveIf(
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

BOOST_AUTO_TEST_CASE(randomTests) {
  BOOST_CHECK(temple::random.getSingle<unsigned>(0, 0) == 0);
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
