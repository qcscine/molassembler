#define BOOST_TEST_MODULE ConstexprMagicTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "template_magic/TemplateMagic.h"
#include "template_magic/Random.h"

#include "Math.h"

#include <iostream>
#include <iomanip>

using namespace std::string_literals;

BOOST_AUTO_TEST_CASE( mathApproxEqual ) {
  const unsigned numTests = 100;

  constexpr double accuracy = 1e-12;

  static_assert(
    accuracy >= std::numeric_limits<double>::epsilon(),
    "Testing accuracy must be greater than machine epsilon!"
  );

  // sqrt
  const auto randomPositiveNumbers = TemplateMagic::random.getN<double>(0, 1e6, numTests);
  auto sqrt_passes = TemplateMagic::map(
    randomPositiveNumbers,
    [](const double& randomPositiveNumber) -> bool {
      return ConstexprMagic::Math::isCloseRelative(
        ConstexprMagic::Math::sqrt(randomPositiveNumber),
        std::sqrt(randomPositiveNumber),
        accuracy
      );
    }
  );

  bool all_sqrt_pass = TemplateMagic::all_of(sqrt_passes);

  BOOST_CHECK(all_sqrt_pass);

  if(!all_sqrt_pass) {
    std::cout << "Square-root implementation is lacking! Failures: " << std::endl;

    for(unsigned i = 0; i < sqrt_passes.size(); i++) {
      if(!sqrt_passes[i]) {
        std::cout << "  x = " << std::setw(12) << randomPositiveNumbers[i]
          << ", sqrt = " << std::setw(12) << ConstexprMagic::Math::sqrt(randomPositiveNumbers[i])
          << ", std::sqrt = " << std::setw(12) << std::sqrt(randomPositiveNumbers[i])
          << ", |Δ| = " << std::setw(12) << std::fabs(
            ConstexprMagic::Math::sqrt(randomPositiveNumbers[i])
            - std::sqrt(randomPositiveNumbers[i])
          ) << std::endl;
      }
    }

    std::cout << std::endl;
  }

  // asin
  const auto randomInverseTrigNumbers = TemplateMagic::random.getN<double>(
    -1 + std::numeric_limits<double>::epsilon(),
    1 - std::numeric_limits<double>::epsilon(),
    numTests
  );

  auto asin_passes = TemplateMagic::map(
    randomInverseTrigNumbers,
    [](const double& randomInverseTrigNumber) -> bool {
      return ConstexprMagic::Math::isCloseRelative(
        ConstexprMagic::Math::asin(randomInverseTrigNumber),
        std::asin(randomInverseTrigNumber),
        accuracy
      );
    }
  );

  bool all_asin_pass = TemplateMagic::all_of(asin_passes);

  BOOST_CHECK(all_asin_pass);

  if(!all_asin_pass) {
    std::cout << "Inverse sine implementation is lacking! Failures: " << std::endl;

    for(unsigned i = 0; i < asin_passes.size(); i++) {
      if(!asin_passes[i]) {
        std::cout << "  x = " << std::setw(12) << randomInverseTrigNumbers[i]
          << ", asin = " << std::setw(12) << ConstexprMagic::Math::asin(randomInverseTrigNumbers[i])
          << ", std::asin = " << std::setw(12) << std::asin(randomInverseTrigNumbers[i])
          << ", |Δ| = " << std::setw(12) << std::fabs(
            ConstexprMagic::Math::asin(randomInverseTrigNumbers[i])
            - std::asin(randomInverseTrigNumbers[i])
          ) << std::endl;
      }
    }

    std::cout << std::endl;
  }

  auto testPow = [](const double& number, const int& exponent) -> bool {
    const double test = ConstexprMagic::Math::pow(number, exponent);
    const double reference = std::pow(number, exponent);

    bool passes = ConstexprMagic::Math::isCloseRelative(
      test,
      reference,
      accuracy
    );

    if(!passes) {
      std::cout << "  x = " << std::setw(12) << number
        << ", exp = " << std::setw(4) << exponent
        << ", pow = " << std::setw(12) << test
        << ", std::pow = " << std::setw(12) << reference
        << ", |Δ| = " << std::setw(12) << std::fabs(test - reference) << ", max permissible diff: "
        << (
          accuracy * std::max(
            std::fabs(test),
            std::fabs(reference)
          )
        ) << std::endl;
    }

    return passes;
  };

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::zipMapAlternate(
        TemplateMagic::random.getN<double>(-1e5, 1e5, numTests),
        TemplateMagic::random.getN<int>(-40, 40, numTests),
        testPow
      )
    )
  );

  // ln
  const auto randomZ = TemplateMagic::random.getN<double>(1e-10, 1e10, numTests);
  bool all_ln_pass = TemplateMagic::all_of(
    TemplateMagic::map(
      randomZ,
      [](const auto& z) -> bool {
        bool pass = ConstexprMagic::Math::isCloseRelative(
          ConstexprMagic::Math::ln(z),
          std::log(z),
          accuracy
        );

        if(!pass) {
          std::cout << "ln deviates for z = " << std::setw(12) << z 
            << ", ln(z) = " << std::setw(12) << ConstexprMagic::Math::ln(z) 
            << ", std::log(z) = " << std::setw(12) << std::log(z)
            << ", |Δ| = " << std::setw(12) << std::fabs(
              ConstexprMagic::Math::ln(z) - std::log(z)
            ) << std::endl;
        }

        return pass;
      }
    )
  );

  BOOST_CHECK(all_ln_pass);

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::map(
        TemplateMagic::random.getN<double>(-100, 100, numTests),
        [](const double& x) -> bool {
          return(ConstexprMagic::Math::floor(x) <= x);
        }
      )
    )
  );

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::map(
        TemplateMagic::random.getN<double>(-100, 100, numTests),
        [](const double& x) -> bool {
          return(ConstexprMagic::Math::ceil(x) >= x);
        }
      )
    )
  );

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::map(
        TemplateMagic::random.getN<double>(-100, 100, numTests),
        [](const double& x) -> bool {
          const double rounded = ConstexprMagic::Math::round(x);
          return(
            rounded == ConstexprMagic::Math::floor(x)
            || rounded == ConstexprMagic::Math::ceil(x)
          );
        }
      )
    )
  );
}
