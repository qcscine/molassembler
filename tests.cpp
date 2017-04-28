#define BOOST_TEST_MODULE ConnectivityManagerTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "templateMagic.h"
#include "Enumerate.h"

// TEMPORARY
#include <iostream>

#include <vector>

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

/*BOOST_AUTO_TEST_CASE( enumerateTests) {
  std::vector<unsigned> testVec {5, 2, 3, 4};

  std::cout << "Before enumerate:" << std::endl;
  for(const auto& enumStruct : enumerate(testVec)) {
    std::cout << "{ index: " << enumStruct.index << ", value: " 
      << enumStruct.value << "}" << std::endl;
  }
}*/
