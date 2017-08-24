#define BOOST_TEST_MODULE ConnectivityManagerTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Enumerate.h"
#include "Random.h"
#include "Containers.h"
#include "Numeric.h"
#include "MemberFetcher.h"
#include "VectorView.h"

// TEMPORARY
#include <iostream>

#include <vector>
#include <cmath>

#include <chrono>

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

  return TemplateMagic::average(timings);
}

double divByThree (unsigned a) {
  return static_cast<double>(a) / 3.0;
}

BOOST_AUTO_TEST_CASE( sumTest ) {
  std::vector<unsigned> instance {0, 1, 2, 3};
  auto f = TemplateMagic::sum(instance);

  BOOST_CHECK(f == 6);

  auto mapped = TemplateMagic::map(
    instance,
    divByThree
  );

  BOOST_CHECK(mapped == std::vector<double>({0, 1.0/3.0, 2.0/3.0, 1}));

  auto pairwiseSum = TemplateMagic::mapSequentialPairs(
    instance,
    std::plus<unsigned>()
  );

  BOOST_CHECK(pairwiseSum == std::vector<unsigned>({1,3,5}));

  auto pairwiseSmaller = TemplateMagic::accumulate(
    TemplateMagic::mapSequentialPairs(
      instance,
      std::less<unsigned>()
    ),
    true,
    std::logical_and<bool>()
  );

  BOOST_CHECK(pairwiseSmaller);

  std::vector<
    std::vector<unsigned>
  > vectorOfVectors {
    {0, 1, 4},
    {4, 5}
  };

  auto mapToSizes = TemplateMagic::map(
    vectorOfVectors,
    [](const std::vector<unsigned>& vectorUnsigned) -> unsigned {
      return vectorUnsigned.size();
    }
  );

  std::vector<unsigned> unsignedVector {1, 2, 3};

  BOOST_CHECK(
    TemplateMagic::sum(
      TemplateMagic::mapAllPairs(
        unsignedVector,
        [](const unsigned& a, const unsigned& b) -> unsigned {
          return a + b;
        }
      )
    ) == 12
  );

  std::vector<double> doubleVector {1.2, 1.5, 1.9};

  BOOST_CHECK(
    TemplateMagic::sum(
      TemplateMagic::mapAllPairs(
        doubleVector,
        [](const double& a, const double& b) -> double {
          return a + b;
        }
      )
    ) == 9.2
  );
}

BOOST_AUTO_TEST_CASE( reduceTests) {
  std::vector<unsigned> values {1, 2, 3, 4, 5};
  BOOST_CHECK(
    TemplateMagic::reduce(
      values,
      0u,
      std::plus<unsigned>()
    ) == 15u
  );
  BOOST_CHECK(
    TemplateMagic::reduce(
      values,
      1u,
      std::multiplies<unsigned>()
    ) == 120u
  );
}

BOOST_AUTO_TEST_CASE( minMaxTests ) {
  const std::vector<unsigned> values {1, 4, 6, 8};
  BOOST_CHECK(TemplateMagic::max(values) == 8u);
  BOOST_CHECK(TemplateMagic::min(values) == 1u);
}

BOOST_AUTO_TEST_CASE(kahanSummation) {
  for(unsigned nTest = 0; nTest < 100; nTest++) {
    const unsigned N = 100;
    const unsigned magnitudeSpread = 20;
    const auto randomNumbers = TemplateMagic::random.getN<double>(
      std::pow(10, - static_cast<double>(magnitudeSpread) / 2),
      std::pow(10, static_cast<double>(magnitudeSpread) / 2),
      N
    );

    const double reduceSum = TemplateMagic::reduce(
      randomNumbers,
      0.0,
      std::plus<double>()
    );

    const double kahanSum = TemplateMagic::kahanSum(randomNumbers);

    /* Reference sum with long doubles, I know an alternative implementation of
     * standard reduce summation with long intermediates would have done the
     * trick too but I'm lazy and this is just a test
     */
    const auto addedPrecision = TemplateMagic::cast<long double>(randomNumbers);
    const double longSum = TemplateMagic::reduce(
      addedPrecision,
      0.0l,
      std::plus<long double>()
    );

    // Kahan summation should be equally or more accurate than the standard reduction
    BOOST_CHECK_MESSAGE(
      std::fabs(longSum - reduceSum) >= std::fabs(longSum - kahanSum),
      "Kahan summation is less accurate than standard reduce sum! "
      << "long: " << longSum << ", kahan: " << kahanSum 
      << ", reduce: " << reduceSum << ". Absolute deviation from long sum: "
      << "kahan:" << std::fabs(longSum - kahanSum) << ", "
      << "reduce: " << std::fabs(longSum - reduceSum)
    );
  }
}

BOOST_AUTO_TEST_CASE(numericAverageStdDev) {
  const std::vector<double> values {29, 30, 31, 32, 33};

  BOOST_CHECK(TemplateMagic::average(values) == 31); 
  BOOST_CHECK(
    std::fabs(TemplateMagic::stddev(values) - std::sqrt(2)) 
    < 1e-10
  );
}

BOOST_AUTO_TEST_CASE(mapToSameContainerTests) {
  std::set<int> f {5, -1, 9};

  auto fMapped = TemplateMagic::map(
    f,
    [](const int& x) -> double {
      return x + 1.3;
    }
  );

  static_assert(
    std::is_same<decltype(fMapped), std::set<double>>::value, 
    "Map to same container does not work as expected"
  );

  std::vector<float> x {0, 3.4, 9};
  auto xMapped = TemplateMagic::map(
    x,
    [](const float& x) -> unsigned long {
      return static_cast<unsigned long>(x);
    }
  );

  static_assert(
    std::is_same<decltype(xMapped), std::vector<unsigned long>>::value, 
    "Map to same container does not work as expected"
  );
}

BOOST_AUTO_TEST_CASE( boolArrayAllOf) {
  std::array<bool, 3> testArray {true, false, true};

  BOOST_CHECK(!TemplateMagic::all_of(testArray));
}

/*BOOST_AUTO_TEST_CASE( enumerateTests) {
  std::vector<unsigned> testVec {5, 2, 3, 4};

  std::cout << "Before enumerate:" << std::endl;
  for(const auto& enumStruct : enumerate(testVec)) {
    std::cout << "{ index: " << enumStruct.index << ", value: " 
      << enumStruct.value << "}" << std::endl;
  }
}*/

BOOST_AUTO_TEST_CASE(memberFetcherTests) {
  std::vector<
    std::vector<unsigned>
  > testClass {
    {4, 3, 9, 10},
    {1, 2, 3, 4, 5},
    {11},
    {44, 93}
  };

  auto memberFetcherClass = TemplateMagic::getMember(
    testClass,
    [](const auto& container) -> unsigned {
      return container.size();
    }
  );

  BOOST_CHECK_MESSAGE(
    TemplateMagic::sum(memberFetcherClass) == 12u,
    "Member fetcher does not work as expected, begin()-end() contains {"
      << TemplateMagic::condenseIterable(memberFetcherClass)
      << "} instead of expected {4, 5, 1, 2}"
  );

  BOOST_CHECK_MESSAGE(
    TemplateMagic::min(memberFetcherClass) == 1
    && TemplateMagic::max(memberFetcherClass) == 5,
    "Usage of min/max does not return expected result"
  );

  // Is this approach more efficient than copying member information?
  double copyingTime = timeNullaryCallable(
    [&]() {
      auto mappedData = TemplateMagic::map(
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
      auto fetcher = TemplateMagic::getMember(
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

  auto filteredView = TemplateMagic::filter(
    testClass,
    [](const auto& subVector) -> bool {
      return subVector.size() < 3;
    }
  );

  auto minSize = TemplateMagic::min(
    TemplateMagic::getMember(
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

  auto maxSize = TemplateMagic::max(
    TemplateMagic::getMember(
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
