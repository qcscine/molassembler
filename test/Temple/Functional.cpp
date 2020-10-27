/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Adaptors/All.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/constexpr/Numeric.h"

// boost algorithm replacements
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/algorithm/transform.hpp>

using namespace Scine::Molassembler;

PURITY_STRONG double divByThree (unsigned a) {
  return static_cast<double>(a) / 3.0;
}

BOOST_AUTO_TEST_CASE(SumTest, *boost::unit_test::label("Temple")) {
  std::vector<unsigned> instance {0, 1, 2, 3};
  auto f = Temple::sum(instance);

  BOOST_CHECK(f == 6);

  BOOST_CHECK(
    boost::accumulate(instance, 0) == 6
  );

  auto mapped = Temple::map(instance, divByThree);
  BOOST_CHECK(mapped == std::vector<double>({0, 1.0/3.0, 2.0/3.0, 1}));

  std::vector<double> bmapped;
  boost::transform(instance, std::back_inserter(bmapped), divByThree);
  BOOST_CHECK(mapped == bmapped);

  auto pairwiseSum = Temple::map(
    Temple::Adaptors::sequentialPairs(instance),
    std::plus<>()
  );

  BOOST_CHECK(pairwiseSum == std::vector<unsigned>({1,3,5}));

  std::vector<unsigned> bpairwiseSum;
  boost::transform(
    instance,
    boost::adaptors::slice(instance, 1, instance.size() - 1),
    std::back_inserter(bpairwiseSum),
    std::plus<>()
  );

  auto pairwiseSmaller = Temple::accumulate(
    Temple::map(
      Temple::Adaptors::sequentialPairs(instance),
      std::less<>()
    ),
    true,
    std::logical_and<>()
  );

  auto bpairwiseSmaller = boost::algorithm::all_of(
    boost::combine(
      boost::adaptors::slice(instance, 0, instance.size() - 2),
      boost::adaptors::slice(instance, 1, instance.size() - 1)
    ),
    [&](const auto& twoTuple) -> bool {
      return Temple::invoke(std::less<>(), twoTuple);
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

  auto mapToSizes = Temple::map(
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

  BOOST_CHECK(
    Temple::sum(
      Temple::map(
        Temple::Adaptors::allPairs(unsignedVector),
        [](const unsigned a, const unsigned b) -> unsigned {
          return a + b;
        }
      )
    ) == 12
  );

  BOOST_CHECK(
    !Temple::all_of(
      unsignedVector,
      [](auto i) -> bool {
        return i < 2;
      }
    )
  );

  std::vector<double> doubleVector {1.2, 1.5, 1.9};

  BOOST_CHECK(
    Temple::sum(
      Temple::map(
        Temple::Adaptors::allPairs(doubleVector),
        [](const double a, const double b) -> double {
          return a + b;
        }
      )
    ) == 9.2
  );
}

BOOST_AUTO_TEST_CASE(ReduceTests, *boost::unit_test::label("Temple")) {
  std::vector<unsigned> values {1, 2, 3, 4, 5};
  BOOST_CHECK(
    Temple::accumulate(
      values,
      0U,
      std::plus<>()
    ) == 15U
  );
  BOOST_CHECK(
    Temple::accumulate(
      values,
      1U,
      std::multiplies<>()
    ) == 120U
  );
}

BOOST_AUTO_TEST_CASE(MinMaxTests, *boost::unit_test::label("Temple")) {
  const std::vector<unsigned> values {1, 4, 6, 8};
  BOOST_CHECK(Temple::max(values) == 8);
  BOOST_CHECK(Temple::min(values) == 1);
}


BOOST_AUTO_TEST_CASE(MapToSameContainerTests, *boost::unit_test::label("Temple")) {
  std::set<int> f {5, -1, 9};

  auto fMapped = Temple::map_stl(
    f,
    [](const int& x) -> double {
      return x + 1.3;
    }
  );

  static_assert(
    std::is_same<decltype(fMapped), std::set<double>>::value,
    "Map to same container does not work as expected"
  );

  std::vector<float> values {0, 3.4, 9};
  auto xMapped = Temple::map(
    values,
    [](const float& x) -> unsigned long {
      return static_cast<unsigned long>(x);
    }
  );

  static_assert(
    std::is_same<decltype(xMapped), std::vector<unsigned long>>::value,
    "Map to same container does not work as expected"
  );
}

BOOST_AUTO_TEST_CASE(MapTuplelikes, *boost::unit_test::label("Temple")) {
  auto f = Temple::map(std::make_pair(4, -1), [](const int x) { return x > 0; });
  static_assert(std::is_same<decltype(f), std::pair<bool, bool>>::value, "Different return type than expected");

  auto g = Temple::map(std::make_tuple(4U, -9, 5.9), [](auto x) { return x + 1; });
  static_assert(std::is_same<decltype(g), std::tuple<unsigned, int, double>>::value, "different return type than expected");

  auto h = Temple::map(std::array<int, 2> {{-1, 4}}, [](int x) { return x > 0; });
  static_assert(std::is_same<decltype(h), std::array<bool, 2>>::value, "Different return type than expected");
}
