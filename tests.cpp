#define BOOST_TEST_MODULE CyclicPolygonsTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "CyclicPolygons.h"
#include <iostream>
#include "template_magic/VectorView.h"
#include "template_magic/Random.h"

BOOST_AUTO_TEST_CASE(symmetricPolynomialsCorrect) {
  std::vector<double> values {1, 2, 3, 4, 5};
  std::vector<double> ks {0, 1, 2, 3, 4, 5};
  std::vector<double> expectedResults {1, 15, 85, 225, 274, 120};

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::zipMapAlternate(
        ks,
        expectedResults,
        [&values](const double& k, const double& expectedValue) -> bool {
          return (
            CyclicPolygons::math::elementarySymmetricPolynomial(
              k,
              values
            ) == expectedValue
          );
        }
      )
    )
  );
}

BOOST_AUTO_TEST_CASE(intSeqTests) {
  std::vector<unsigned> sequence {1, 2, 3, 4, 5};

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::zipMapAlternate(
        sequence,
        CyclicPolygons::math::intSeq(sequence.front(), sequence.back()),
        std::equal_to<unsigned>()
      )
    )
  );
}

BOOST_AUTO_TEST_CASE(comparisonAgainstRImplementation) {
  std::vector<double> edgeLengths {29, 30, 31, 32, 33};

  const auto squaredEdgeLengths = TemplateMagic::map(
    edgeLengths,
    CyclicPolygons::math::square<double>
  );

  const auto epsilon = TemplateMagic::map(
    CyclicPolygons::math::intSeq(0, 5),
    [&](const unsigned& k) -> double {
      return CyclicPolygons::math::elementarySymmetricPolynomial(
        k,
        squaredEdgeLengths
      );
    }
  );
  
  // Compare epsilons
  const std::vector<double> rEpsilons {
    1.000000e+00,
    4.815000e+03,
    9.254463e+06,
    8.875070e+09,
    4.246737e+12,
    8.111286e+14
  };

  auto relativeEquals = [](const double& relativeTolerance) {
    return [&](const double& calc, const double& ref) -> bool {
      if(std::fabs(calc - ref) > relativeTolerance * std::fabs(std::max(calc, ref))) {
        std::cout << "Values do not match to within " << relativeTolerance
          << ": " << calc << ", reference: " << ref << std::endl;
        /*std::cout << "fabs diff:" << std::fabs(calc - ref) << " > " 
          << relativeTolerance << " * " << std::fabs(std::max(calc, ref)) 
          << " (fabs-max)" << std::endl;*/

        return false;
      }

      return true;
    };
  };

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::zipMapAlternate(
        epsilon,
        rEpsilons,
        relativeEquals(1e-6)
      )
    )
  );

  // Lambda for a given rho
  const double rho = 0.01;
  const auto lambdas = TemplateMagic::map(
    CyclicPolygons::math::intSeq(0, 4),
    [&](const unsigned& k) -> double {
      return CyclicPolygons::Pentagon::lambda(
        k,
        rho,
        epsilon
      );
    }
  );

  const std::vector<double> rLambdas {
    -34038.1109,
    18362.3869, 
    -4550.5927,
    585.2463,
    -38.1500  
  };

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::zipMapAlternate(
        lambdas,
        rLambdas,
        relativeEquals(1e-6)
      )
    )
  );

  const double rA5 = 1044227;
  const double rB5 = -90233.8;
  const double rDelta5 = -5318.33;

  BOOST_CHECK(
    relativeEquals(1e-6)(
      CyclicPolygons::Pentagon::A5(lambdas),
      rA5
    )
  );
  BOOST_CHECK(
    relativeEquals(1e-6)(
      CyclicPolygons::Pentagon::B5(lambdas),
      rB5
    )
  );
  BOOST_CHECK(
    relativeEquals(1e-6)(
      CyclicPolygons::Pentagon::Delta5(lambdas),
      rDelta5
    )
  );
}

BOOST_AUTO_TEST_CASE(findsCorrectRoot) {
  using namespace std::string_literals;

  /* Limited case, where shortest and longest possible bond lengths in molecules
   * are used
   */
  const double upperLimit = 5.6; // Fr-Fr single
  const double lowerLimit = 1.6; // H-H single

  const unsigned nTests = 100;

  unsigned nFailures = 0;
  for(unsigned nTest = 0; nTest < nTests; ++nTest) {
    std::vector<double> edgeLengths {
      TemplateMagic::random.getSingle<double>(
        lowerLimit,
        upperLimit
      )
    };

    while(edgeLengths.size() < 5) {
      const double geometricLimit = TemplateMagic::numeric::sum(edgeLengths);

      edgeLengths.emplace_back(
        TemplateMagic::random.getSingle<double>(
          lowerLimit,
          std::min(upperLimit, geometricLimit)
        )
      );
    }

    auto radiusOption = CyclicPolygons::Pentagon::maximumCircumradius(edgeLengths);

    // Finds a radius
    if(!radiusOption) {
      nFailures += 1;
      CyclicPolygons::analysis::writeAnalysisFiles(
        edgeLengths,
        "test-failure-"s + std::to_string(nTest)
      );
      std::cout << "For this set of edgeLengths, no valid root was found" << std::endl;
    }

    if(radiusOption) {
      // radius is correct
      const double rho = 1 / std::pow(radiusOption.value(), 2);
      if(!CyclicPolygons::Pentagon::validateRhoGuess(edgeLengths, rho)) {
        nFailures += 1;
        CyclicPolygons::analysis::writeAnalysisFiles(
          edgeLengths,
          "test-failure-"s + std::to_string(nTest)
        );
        std::cout << "maximumPentagonCircumradius returned non-validated "
          << "circumradius!" << std::endl;
      }
    }
  }

  BOOST_CHECK_MESSAGE(
    nFailures == 0,
    "Out of " << nTests << " root-finding tests, " << nFailures << " failed."
  );
}

BOOST_AUTO_TEST_CASE(internalAnglesSumCorrectly) {
  std::vector<double> edgeLengths {29, 30, 31, 32, 33};
  for(unsigned n = 3; n <= 5; n++) {
    auto edges = TemplateMagic::cast<double>(
      CyclicPolygons::math::intSeq(29, 28 + n)
    );

    auto internalAngles = CyclicPolygons::internalAngles(edges);

    BOOST_CHECK(
      std::fabs(
        TemplateMagic::numeric::sum(internalAngles)
        - M_PI * (n - 2)
      ) < 1e-6
    );
  }
}
