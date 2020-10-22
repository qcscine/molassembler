/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/results_collector.hpp>
#include <boost/test/unit_test.hpp>

#include "boost/optional.hpp"

#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Functor.h"
#include "Molassembler/Temple/OperatorSuppliers.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/constexpr/Numeric.h"
#include "Molassembler/Temple/constexpr/Optional.h"

#include <cmath>
#include <vector>
#include <array>

using namespace Scine::Molassembler;

template<typename T>
std::ostream& operator << (std::ostream& os, const boost::optional<T>& valueOptional) {
  if(valueOptional) {
    os << "Some " << valueOptional.value();
  } else {
    os << "None";
  }

  return os;
}

BOOST_AUTO_TEST_CASE(SetRemovalDefect, *boost::unit_test::label("Temple")) {
  /* Sets and maps cannot use std::remove_if! Not a defect. */

  auto testSet = Temple::copy_if(
    std::set<unsigned> {5, 2, 9, 1, 3, 4, 8},
    [](const auto& value) -> bool {
      return value % 2 != 0;
    }
  );

  BOOST_CHECK((testSet == std::set<unsigned> {1, 3, 5, 9}));
}

BOOST_AUTO_TEST_CASE(SelectTestCases, *boost::unit_test::label("Temple")) {
  auto ragged2D = std::vector<
    std::vector<unsigned>
  > {
    {4, 9, 3},
    {2, 11, 6},
    {30, 4},
    {9, 44, 33, 12}
  };

  auto smallestVector = Temple::select(
    ragged2D,
    std::less<>(),
    [](const auto& group) -> unsigned {
      return group.size();
    }
  );

  BOOST_CHECK((*smallestVector == std::vector<unsigned> {30, 4}));

  auto largestVector = Temple::select(
    ragged2D,
    std::greater<>(),
    [](const auto& group) -> unsigned {
      return group.size();
    }
  );

  BOOST_CHECK((*largestVector == std::vector<unsigned> {9, 44, 33, 12}));
}

BOOST_AUTO_TEST_CASE(StringifyTests, *boost::unit_test::label("Temple")) {
  std::vector<
    std::map<
      unsigned,
      std::pair<int, double>
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

  std::string stringified = Temple::stringify(complicatedStructure);

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

  stringified = Temple::stringify(someTuple);

  BOOST_CHECK(!stringified.empty());
}

struct TupleLike : Temple::Crtp::LexicographicComparable<TupleLike> {
  unsigned x;
  int f;

  TupleLike(unsigned x_, int f_) : x(x_), f(f_) {}

  auto tie() const {
    return std::tie(x, f);
  }
};

BOOST_AUTO_TEST_CASE(CrtpTests, *boost::unit_test::label("Temple")) {
  TupleLike a {4, -3};
  TupleLike b {9, 4};
  TupleLike c {9, 4};
  BOOST_CHECK(b == c);
  BOOST_CHECK(a < b);
  BOOST_CHECK(b > a);
  BOOST_CHECK(b <= c);
  BOOST_CHECK(b >= c);
}

BOOST_AUTO_TEST_CASE(ReferenceOptional, *boost::unit_test::label("Temple")) {
  int x = 0;
  Temple::Optional<int&> f {x};
  if(f) {
    f.value() = 4;
  }
  BOOST_CHECK(x == 4);

  Temple::Optional<int&> g;
  BOOST_CHECK(!g.hasValue());
}

BOOST_AUTO_TEST_CASE(FunctorSafety, *boost::unit_test::label("Temple")) {
  std::vector<unsigned> x {1, 4, 3};
  // Bind x by const reference within the functor
  auto at = Temple::Functor::at(x);
  BOOST_CHECK(at(0) == x.at(0));
  // External modification affects the functor too since access is referential
  x.push_back(5);
  BOOST_CHECK(at(3) == x.at(3));

  auto at2 = Temple::Functor::at(std::vector<unsigned> {1, 2, 3});
  BOOST_CHECK(at2(0) == 1);

  constexpr std::tuple<int, double> t {4, 1.0};
  static_assert(Temple::Functor::get<0>()(t) == 4, "Get functor doesn't work");
  static_assert(Temple::Functor::first(t) == 4, "Pair_first doesn't work");
  BOOST_CHECK_EQUAL(Temple::Functor::second(t), 1.0);
}
