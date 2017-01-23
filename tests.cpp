#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include "templateMagic.h"

#include <vector>

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
}
