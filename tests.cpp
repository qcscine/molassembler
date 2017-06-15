#define BOOST_TEST_MODULE ConnectivityManagerTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "TemplateMagic.h"
#include "Enumerate.h"
#include "Random.h"

// TEMPORARY
#include <iostream>

#include <vector>
#include <cmath>

double divByThree (unsigned a) {
  return static_cast<double>(a) / 3.0;
}

BOOST_AUTO_TEST_CASE( sumTest ) {
  std::vector<unsigned> instance {0, 1, 2, 3};
  auto f = TemplateMagic::numeric::sum(instance);

  BOOST_CHECK(f == 6);

  auto mapped = TemplateMagic::map(
    instance,
    divByThree
  );

  BOOST_CHECK(mapped == std::vector<double>({0, 1.0/3.0, 2.0/3.0, 1}));

  auto pairwiseSum = TemplateMagic::pairwiseMap(
    instance,
    std::plus<unsigned>()
  );

  BOOST_CHECK(pairwiseSum == std::vector<unsigned>({1,3,5}));

  auto pairwiseSmaller = TemplateMagic::accumulate(
    TemplateMagic::pairwiseMap(
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
    TemplateMagic::numeric::sum(
      TemplateMagic::allPairsMap(
        unsignedVector,
        [](const unsigned& a, const unsigned& b) -> unsigned {
          return a + b;
        }
      )
    ) == 12
  );

  std::vector<double> doubleVector {1.2, 1.5, 1.9};

  BOOST_CHECK(
    TemplateMagic::numeric::sum(
      TemplateMagic::allPairsMap(
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
  BOOST_CHECK(TemplateMagic::numeric::max(values) == 8u);
  BOOST_CHECK(TemplateMagic::numeric::min(values) == 1u);
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

    const double kahanSum = TemplateMagic::numeric::kahanSum(randomNumbers);

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

  BOOST_CHECK(TemplateMagic::numeric::average(values) == 31); 
  BOOST_CHECK(
    std::fabs(TemplateMagic::numeric::stddev(values) - std::sqrt(2)) 
    < 1e-10
  );
}

/*BOOST_AUTO_TEST_CASE( enumerateTests) {
  std::vector<unsigned> testVec {5, 2, 3, 4};

  std::cout << "Before enumerate:" << std::endl;
  for(const auto& enumStruct : enumerate(testVec)) {
    std::cout << "{ index: " << enumStruct.index << ", value: " 
      << enumStruct.value << "}" << std::endl;
  }
}*/
