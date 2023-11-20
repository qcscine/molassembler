/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/constexpr/Math.h"
#include "Molassembler/Temple/constexpr/Jsf.h"
#include "Molassembler/Temple/constexpr/FloatingPointComparison.h"
#include "Molassembler/Temple/Stringify.h"

#include <iomanip>
#include <iostream>

using namespace Scine::Molassembler;
extern Temple::Generator<> generator;

static_assert(Temple::Math::factorial(5) == 120, "Factorial is incorrect");
static_assert(Temple::Math::factorial(0) == 1, "Factorial is incorrect");

namespace {

constexpr unsigned numTests = 100;
constexpr double relativeAccuracy = 1e-12;

static_assert(
  relativeAccuracy >= std::numeric_limits<double>::epsilon(),
  "Testing relative accuracy must be greater than machine epsilon!"
);

template<typename F, typename G>
auto compareImplFn(F&& testFn, G&& referenceFn, double accuracy=relativeAccuracy) {
  return [=](const auto ... values) {
    const double testValue = Temple::invoke(testFn, values...);
    const double referenceValue = Temple::invoke(referenceFn, values...);

    BOOST_TEST_CONTEXT(
      "  x = " << std::setw(12) << Temple::stringify(std::tie(values...))
      << ", y = " << std::setw(12) << testValue
      << ", ref = " << std::setw(12) << referenceValue
      << ", |Δ| = " << std::setw(12) << std::fabs(testValue - referenceValue)
    ) {
      BOOST_CHECK(
        Temple::Floating::isCloseRelative(testValue, referenceValue, accuracy)
      );
    }
  };
}

} // namespace

BOOST_AUTO_TEST_CASE(ConstexprSqrt, *boost::unit_test::label("Temple")) {
  Temple::forEach(
    Temple::Random::getN<double>(0, 1e6, numTests, generator.engine),
    compareImplFn(
      [](double x) { return Temple::Math::sqrt(x); },
      [](double x) { return std::sqrt(x); }
    )
  );
}

BOOST_AUTO_TEST_CASE(ConstexprAsin, *boost::unit_test::label("Temple")) {
  // asin
  const auto randomInverseTrigNumbers = Temple::Random::getN<double>(
    std::nexttoward(-1.0, 0.0),
    std::nexttoward(1.0, 0.0),
    numTests,
    generator.engine
  );

  Temple::forEach(
    randomInverseTrigNumbers,
    compareImplFn(
      [](double x) { return Temple::Math::asin(x); },
      [](double x) { return std::asin(x); },
      1e-8
    )
  );
}

BOOST_AUTO_TEST_CASE(ConstexprPow, *boost::unit_test::label("Temple")) {
  Temple::forEach(
    Temple::Adaptors::zip(
      Temple::Random::getN<double>(-1e5, 1e5, numTests, generator.engine),
      Temple::Random::getN<int>(-40, 40, numTests, generator.engine)
    ),
    compareImplFn(
      [](double x, int y) { return Temple::Math::pow(x, y); },
      [](double x, int y) { return std::pow(x, y); }
    )
  );
}

BOOST_AUTO_TEST_CASE(ConstexprRecPow, *boost::unit_test::label("Temple")) {
  Temple::forEach(
    Temple::Adaptors::zip(
      Temple::Random::getN<double>(-1e5, 1e5, numTests, generator.engine),
      Temple::Random::getN<unsigned>(0, 40, numTests, generator.engine)
    ),
    compareImplFn(
      [](double x, unsigned y) { return Temple::Math::recPow(x, y); },
      [](double x, unsigned y) { return std::pow(x, y); }
    )
  );
}

BOOST_AUTO_TEST_CASE(ConstexprLn, *boost::unit_test::label("Temple")) {
  Temple::forEach(
    Temple::Random::getN<double>(1e-10, 1e10, numTests, generator.engine),
    compareImplFn(
      [](double x) { return Temple::Math::ln(x); },
      [](double x) { return std::log(x); }
    )
  );
}

BOOST_AUTO_TEST_CASE(ConstexprAtan, *boost::unit_test::label("Temple")) {
  Temple::forEach(
    Temple::Random::getN<double>(-M_PI / 2, M_PI / 2, numTests, generator.engine),
    compareImplFn(
      [](double x) { return Temple::Math::atan(x); },
      [](double x) { return std::atan(x); }
    )
  );
}

BOOST_AUTO_TEST_CASE(ConstexprFloorCeil, *boost::unit_test::label("Temple")) {
  BOOST_CHECK(
    Temple::all_of(
      Temple::Random::getN<double>(-100, 100, numTests, generator.engine),
      [](const double x) -> bool {
        return(Temple::Math::floor(x) <= x);
      }
    )
  );

  BOOST_CHECK(
    Temple::all_of(
      Temple::Random::getN<double>(-100, 100, numTests, generator.engine),
      [](const double x) -> bool {
        return(Temple::Math::ceil(x) >= x);
      }
    )
  );
}

namespace {

template<
  template<typename> class Comparator,
  typename T
> constexpr bool testComparison(const T a, const T b, const T tolerance) {
  Comparator<T> comparator { tolerance };

  return (
    Temple::Math::XOR(
      (
        comparator.isLessThan(a, b)
        && comparator.isMoreThan(b, a)
        && comparator.isUnequal(a, b)
      ),
      (
        comparator.isLessThan(b, a)
        && comparator.isMoreThan(a, b)
        && comparator.isUnequal(a, b)
      ),
      (
        !comparator.isLessThan(a, b)
        && !comparator.isMoreThan(a, b)
        && comparator.isEqual(a, b)
      )
    ) && Temple::Math::XOR(
      comparator.isEqual(a, b),
      comparator.isUnequal(a, b)
    )
  );
}

using namespace Temple::Floating;

static_assert(
  testComparison<ExpandedAbsoluteEqualityComparator>(4.3, 3.9, 1e-4)
  && testComparison<ExpandedAbsoluteEqualityComparator>(4.3, 3.9, 1.0)
  && testComparison<ExpandedAbsoluteEqualityComparator>(4.4, 4.4, 1e-10),
  "absolute comparison has inconsistent operators!"
);

static_assert(
  testComparison<ExpandedRelativeEqualityComparator>(4.3, 3.9, 1e-4)
  && testComparison<ExpandedRelativeEqualityComparator>(4.3, 3.9, 1.0)
  && testComparison<ExpandedRelativeEqualityComparator>(4.4, 4.4, 1e-10),
  "relative comparison has inconsistent operators!"
);

} // namespace
