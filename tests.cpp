#define BOOST_TEST_MODULE CyclicPolygonsTests
#define BOOST_TEST_DYN_LINK

#include "boost/test/unit_test.hpp"
#include "temple/Containers.h"
#include "temple/Random.h"
#include "temple/VectorView.h"

#include "CyclicPolygons.h"

#include <iostream>

BOOST_AUTO_TEST_CASE(symmetricPolynomialsCorrect) {
  std::vector<double> values {1, 2, 3, 4, 5};
  std::vector<double> ks {0, 1, 2, 3, 4, 5};
  std::vector<double> expectedResults {1, 15, 85, 225, 274, 120};

  BOOST_CHECK(
    temple::all_of(
      temple::zipMap(
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
    temple::all_of(
      temple::zipMap(
        sequence,
        CyclicPolygons::math::intSeq(sequence.front(), sequence.back()),
        std::equal_to<unsigned>()
      )
    )
  );
}

BOOST_AUTO_TEST_CASE(comparisonAgainstRImplementation) {
  std::vector<double> edgeLengths {29, 30, 31, 32, 33};

  const auto squaredEdgeLengths = temple::map(
    edgeLengths,
    CyclicPolygons::math::square<double>
  );

  const auto epsilon = temple::map(
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
    temple::all_of(
      temple::zipMap(
        epsilon,
        rEpsilons,
        relativeEquals(1e-6)
      )
    )
  );

  // Lambda for a given rho
  const double rho = 0.01;
  const auto lambdas = temple::map(
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
    temple::all_of(
      temple::zipMap(
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
      temple::random.getSingle<double>(
        lowerLimit,
        upperLimit
      )
    };

    while(edgeLengths.size() < 5) {
      const double geometricLimit = temple::sum(edgeLengths);

      edgeLengths.emplace_back(
        temple::random.getSingle<double>(
          lowerLimit,
          std::min(upperLimit, geometricLimit)
        )
      );
    }

    auto radiusOption = CyclicPolygons::Pentagon::convexCircumradiusSvrtan(edgeLengths);

    // Finds a radius
    if(!radiusOption) {
      nFailures += 1;
      CyclicPolygons::analysis::writeSvrtanAnalysisFiles(
        edgeLengths,
        "svrtan-test-failure-"s + std::to_string(nTest)
      );
      std::cout << "For this set of edgeLengths, no valid root was found" << std::endl;
    }

    if(radiusOption) {
      // radius is correct
      const double rho = 1 / std::pow(radiusOption.value(), 2);
      if(!CyclicPolygons::Pentagon::validateRhoGuess(edgeLengths, rho)) {
        nFailures += 1;
        CyclicPolygons::analysis::writeSvrtanAnalysisFiles(
          edgeLengths,
          "svrtan-test-failure-"s + std::to_string(nTest)
        );
        std::cout << "maximumPentagonCircumradius returned non-validated "
          << "circumradius!" << std::endl;
      }
    }
  }

  BOOST_CHECK_MESSAGE(
    nFailures == 0,
    "Out of " << nTests << " heuristic root-finding tests, " << nFailures << " failed."
  );
}

BOOST_AUTO_TEST_CASE(centralAngleRootFinding) {
  const double upperLimit = 5.6; // Fr-Fr single
  const double lowerLimit = 0.7; // H-H single

  const double mean = (upperLimit - lowerLimit) / 2;
  const double lowerStddev = 0.01 * mean;
  const double upperStddev = 4 * mean;
  const unsigned nStddevSteps = 40;
  const unsigned nSamples = 1000;

  const double stepLength = (upperStddev - lowerStddev) / nStddevSteps;
  for(unsigned stepNumber = 0; stepNumber <= nStddevSteps; stepNumber++) {
    double currentStddev = lowerStddev + stepNumber * stepLength;

    std::normal_distribution<double> normalDistribution(mean, currentStddev);

    auto sampleTruncatedNormal = [&]() -> double {
      double sample;
      while(true) {
        sample = normalDistribution(
          temple::random.randomEngine
        );
        if(lowerLimit <= sample && sample <= upperLimit) {
          return sample;
        }
      }
    };

    for(unsigned sampleNumber = 0; sampleNumber < nSamples; sampleNumber++) {
      std::vector<double> edgeLengths;
      while(edgeLengths.size() < 5) {
        edgeLengths.emplace_back(sampleTruncatedNormal());
      }

      while(!CyclicPolygons::exists(edgeLengths)) {
        edgeLengths.clear();
        while(edgeLengths.size() < 5) {
          edgeLengths.emplace_back(sampleTruncatedNormal());
        }
      }

      double circumradius;
      auto assignCircumradius = [&]() -> void {
        circumradius = CyclicPolygons::Pentagon::convexCircumradius(
          edgeLengths
        );
      };

      BOOST_CHECK_NO_THROW(assignCircumradius());
    }
  }
}

BOOST_AUTO_TEST_CASE(internalAnglesSumCorrectly) {
  std::vector<double> edgeLengths {29, 30, 31, 32, 33};
  for(unsigned n = 3; n <= 5; n++) {
    auto edges = temple::cast<double>(
      CyclicPolygons::math::intSeq(29, 28 + n)
    );

    auto internalAngles = CyclicPolygons::internalAngles(edges);

    BOOST_CHECK(
      std::fabs(
        temple::sum(internalAngles)
        - M_PI * (n - 2)
      ) < 1e-6
    );
  }
}
