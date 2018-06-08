#define BOOST_TEST_MODULE ContainersTestsModule
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"
#define INCLUDE_INVOKE_WITH_BOOST_TUPLE
#include "temple/Invoke.h"

// boost algorithm replacements
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/algorithm/transform.hpp>

double divByThree (unsigned a) PURITY_STRONG;
double divByThree (unsigned a) {
  return static_cast<double>(a) / 3.0;
}

BOOST_AUTO_TEST_CASE( sumTest ) {
  std::vector<unsigned> instance {0, 1, 2, 3};
  auto f = temple::sum(instance);

  BOOST_CHECK(f == 6);

  BOOST_CHECK(
    boost::accumulate(instance, 0) == 6
  );

  auto mapped = temple::map(
    instance,
    divByThree
  );

  BOOST_CHECK(mapped == std::vector<double>({0, 1.0/3.0, 2.0/3.0, 1}));

  std::vector<double> bmapped;
  boost::transform(instance, std::back_inserter(bmapped), divByThree);
  BOOST_CHECK(mapped == bmapped);

  auto pairwiseSum = temple::mapSequentialPairs(
    instance,
    std::plus<unsigned>()
  );

  BOOST_CHECK(pairwiseSum == std::vector<unsigned>({1,3,5}));

  std::vector<unsigned> bpairwiseSum;
  boost::transform(
    instance,
    boost::adaptors::slice(instance, 1, instance.size() - 1),
    std::back_inserter(bpairwiseSum),
    std::plus<unsigned>()
  );

  auto pairwiseSmaller = temple::accumulate(
    temple::mapSequentialPairs(
      instance,
      std::less<unsigned>()
    ),
    true,
    std::logical_and<bool>()
  );

  auto bpairwiseSmaller = boost::algorithm::all_of(
    boost::combine(
      boost::adaptors::slice(instance, 0, instance.size() - 2),
      boost::adaptors::slice(instance, 1, instance.size() - 1)
    ),
    [&](const auto& twoTuple) -> bool {
      return temple::invoke(std::less<>(), twoTuple);
    }
  );

  static_assert(std::is_same<decltype(bpairwiseSmaller), bool>::value, "Not a bool??");
  BOOST_CHECK(bpairwiseSmaller);

  BOOST_CHECK(pairwiseSmaller);

  std::vector<
    std::vector<unsigned>
  > vectorOfVectors {
    {0, 1, 4},
    {4, 5}
  };

  auto mapToSizes = temple::map(
    vectorOfVectors,
    [](const std::vector<unsigned>& vectorUnsigned) -> unsigned {
      return vectorUnsigned.size();
    }
  );

  std::vector<unsigned> bsizes;
  boost::transform(
    vectorOfVectors,
    std::back_inserter(bsizes),
    [](const std::vector<unsigned>& vectorUnsigned) -> unsigned {
      return vectorUnsigned.size();
    }
  );
  BOOST_CHECK(bsizes == mapToSizes);

  std::vector<unsigned> unsignedVector {1, 2, 3};

  // TODO this one is really difficult to do in boost I think
  BOOST_CHECK(
    temple::sum(
      temple::mapAllPairs(
        unsignedVector,
        [](const unsigned& a, const unsigned& b) -> unsigned {
          return a + b;
        }
      )
    ) == 12
  );

  std::vector<double> doubleVector {1.2, 1.5, 1.9};

  BOOST_CHECK(
    temple::sum(
      temple::mapAllPairs(
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
    temple::reduce(
      values,
      0u,
      std::plus<unsigned>()
    ) == 15u
  );
  BOOST_CHECK(
    temple::reduce(
      values,
      1u,
      std::multiplies<unsigned>()
    ) == 120u
  );
}
BOOST_AUTO_TEST_CASE( minMaxTests ) {
  const std::vector<unsigned> values {1, 4, 6, 8};
  BOOST_CHECK(temple::max(values) == 8u);
  BOOST_CHECK(temple::min(values) == 1u);
}


BOOST_AUTO_TEST_CASE(mapToSameContainerTests) {
  std::set<int> f {5, -1, 9};

  auto fMapped = temple::map(
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
  auto xMapped = temple::map(
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
  std::array<bool, 3> testArray {{true, false, true}};

  BOOST_CHECK(!temple::all_of(testArray));
}

BOOST_AUTO_TEST_CASE(concatenateTests) {
  std::vector<unsigned> f {{4, 9, 1}};
  std::set<unsigned> x {{3, 6, 10}};

  std::vector<unsigned> comp {{4, 9, 1, 3, 6, 10}};

  static_assert(
    std::is_same<
      decltype(temple::concatenate(f, x)),
      std::vector<unsigned>
    >::value,
    "Concatenate must return a vector of the underlying shared type"
  );

  BOOST_CHECK(
    temple::concatenate(f, x) == comp
  );
}

unsigned quaternarySum(
  const unsigned a,
  const unsigned b,
  const unsigned c,
  const unsigned d
) PURITY_STRONG;

unsigned quaternarySum(
  const unsigned a,
  const unsigned b,
  const unsigned c,
  const unsigned d
) {
  return a + b + c + d;
}

BOOST_AUTO_TEST_CASE(unpackArrayTests) {
  std::array<unsigned, 4> values {{1, 2, 3, 5}};
  BOOST_CHECK(
    temple::unpackArrayToFunction(values, quaternarySum) == 11u
  );
}
