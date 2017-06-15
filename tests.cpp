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
      return ConstexprMagic::Math::detail::isApprox(
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
      return ConstexprMagic::Math::detail::isApprox(
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

  // pow
  const auto randomNumbers = TemplateMagic::random.getN<double>(-1e5, 1e5, numTests);
  const auto randomExponents = TemplateMagic::random.getN<int>(-40, 40, numTests);

  std::vector<bool> pow_passes;
  for(unsigned i = 0; i < numTests; i++) {
    pow_passes.emplace_back(
      ConstexprMagic::Math::detail::isApprox(
        ConstexprMagic::Math::pow(randomNumbers[i], randomExponents[i]),
        std::pow(randomNumbers[i], randomExponents[i]),
        accuracy
      )
    );
  }

  bool all_pow_pass = TemplateMagic::all_of(pow_passes);

  if(!all_pow_pass) {
    std::cout << "Power implementation is lacking! Failures: " << std::endl;

    for(unsigned i = 0; i < pow_passes.size(); i++) {
      if(!pow_passes[i]) {
        std::cout << "  x = " << std::setw(12) << randomNumbers[i]
          << ", exp = " << std::setw(4) << randomExponents[i]
          << ", pow = " << std::setw(12) << ConstexprMagic::Math::pow(randomNumbers[i], randomExponents[i])
          << ", std::pow = " << std::setw(12) << std::pow(randomNumbers[i], randomExponents[i])
          << ", |Δ| = " << std::setw(12) << std::fabs(
            ConstexprMagic::Math::pow(randomNumbers[i], randomExponents[i])
            - std::pow(randomNumbers[i], randomExponents[i])
          ) << std::endl;
      }
    }

    std::cout << std::endl;
  }
}
