#define BOOST_TEST_MODULE CyclicPolygonsTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "CyclicPolygons.h"
#include <iostream>
#include <fstream>
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

BOOST_AUTO_TEST_CASE(pentagonCircumradiusTests) {
  std::vector<double> edgeLengths {29, 30, 31, 32, 33};
  std::vector<double> roots {26.385, 16.512, 17.026, 17.595, 17.991, 18.335, 18.651};

  auto rhoRoots = TemplateMagic::map(
    roots,
    [](const double& r) -> double {
      return 1 / (r*r);
    }
  );

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

  std::ofstream outFile("scan.csv");
  outFile << std::scientific << std::setprecision(8);

  const double lower = 0;
  const double upper = 0.0037;
  const unsigned nSteps = 1000;

  const double stepSize = (upper - lower) / nSteps;

  for(unsigned i = 0; i <= nSteps; i++) {
    const double rho = lower + i * stepSize;
    outFile << rho << "," << CyclicPolygons::svrtan::pentagon(
      rho,
      epsilon
    ) << std::endl;
  }

  outFile.close();
}

BOOST_AUTO_TEST_CASE(scanRandomEdgeCombinations) {
  using namespace std::string_literals;

  /* Limited case, where shortest and longest possible bond lengths in molecules
   * are used
   */
  const double upperLimit = 5.6; // Fr-Fr single
  const double lowerLimit = 0.7; // H-H single

  for(unsigned nTest = 0; nTest < 10; ++nTest) {
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

    const auto unaryPolynomial = [&](const double& rho) -> double {
      return CyclicPolygons::svrtan::pentagon(rho, epsilon);
    };

    std::string scanFilename = "scan-random-"s + std::to_string(nTest) + ".csv"s;
    std::string metaFilename = "scan-random-"s + std::to_string(nTest) + "-meta.csv"s;

    std::ofstream outFile(scanFilename);
    outFile << std::scientific << std::setprecision(8);

    const double trialFactor = 1.2;

    const double average = TemplateMagic::numeric::average(edgeLengths);
    const double stddev = TemplateMagic::numeric::stddev(edgeLengths, average);

    const double rhoGuess = 1 / CyclicPolygons::math::square(
      0.5 * average / std::sin(M_PI / 5)
    );
    const double rhoStddev = 2 * stddev / (
      std::pow(rhoGuess, 3) * std::sin(M_PI / edgeLengths.size())
    );
    const bool triedFirstStrategy = (stddev < average * 0.1);

    boost::optional<double> rhoRoot;

    if(triedFirstStrategy) {
      const double lowerValue = unaryPolynomial(rhoGuess / trialFactor);
      const double upperValue = unaryPolynomial(rhoGuess * trialFactor);
      bool increasing = lowerValue < upperValue;

      const unsigned maxIter = 1000;
      boost::uintmax_t iterCount = maxIter; // maxIter
      auto returnBounds = boost::math::tools::bracket_and_solve_root(
        unaryPolynomial,
        rhoGuess,
        trialFactor,
        increasing,
        boost::math::tools::eps_tolerance<double>(32),
        iterCount
      );

      if(iterCount < static_cast<boost::uintmax_t>(maxIter)) {
        // Success! We can return the circumradius
        rhoRoot = (returnBounds.first + returnBounds.second) / 2;
      }
    }

    const bool usedFallbackStrategy = (!rhoRoot);

    const double scanLower = 0;
    const double scanUpper = 2 * rhoGuess;

    const unsigned nSteps = 10000;
    const double stepLength = (scanUpper - scanLower) / nSteps;

    outFile << "0, 0" << std::endl;
    const double initialValue = unaryPolynomial(stepLength);
    outFile << stepLength << ", " << initialValue << std::endl;

    const bool initialSignBit = std::signbit(initialValue);
    for(unsigned i = 2; i < nSteps; i++) {
      const double rho = i * stepLength;
      const double value = unaryPolynomial(rho);

      outFile << rho << ", " << value << std::endl;

      if(
        std::signbit(unaryPolynomial(i * stepLength)) != initialSignBit
        && !rhoRoot
      ) {
        const double lowerBound = (i - 1) * stepLength;
        const double upperBound = rho;

        const double lowerValue = unaryPolynomial(lowerBound);
        const double upperValue = unaryPolynomial(upperBound);
        
        const unsigned maxIter = 1000;
        boost::uintmax_t iterCount = maxIter;

        auto returnBounds = boost::math::tools::toms748_solve(
          unaryPolynomial,
          lowerBound,
          upperBound,
          lowerValue,
          upperValue,
          boost::math::tools::eps_tolerance<double>(32),
          iterCount
        );

        if(iterCount < static_cast<boost::uintmax_t>(maxIter)) {
          // Success! We can return the circumradius
          rhoRoot = (returnBounds.first + returnBounds.second) / 2;
        }
      }
    }

    outFile.close();

    std::ofstream metaFile(metaFilename);
    for(unsigned i = 0; i < edgeLengths.size(); i++) {
      metaFile << edgeLengths[i];
      if(i != edgeLengths.size() - 1) {
        metaFile << ", ";
      }
    }
    metaFile << std::endl;
    metaFile << triedFirstStrategy << ", " << usedFallbackStrategy << std::endl;
    metaFile << rhoGuess << ", " << rhoStddev << std::endl;
    metaFile << rhoRoot.value_or(0) << std::endl;

    metaFile.close();
  }
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
      return CyclicPolygons::svrtan::lambda(
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
      CyclicPolygons::svrtan::A5(lambdas),
      rA5
    )
  );
  BOOST_CHECK(
    relativeEquals(1e-6)(
      CyclicPolygons::svrtan::B5(lambdas),
      rB5
    )
  );
  BOOST_CHECK(
    relativeEquals(1e-6)(
      CyclicPolygons::svrtan::Delta5(lambdas),
      rDelta5
    )
  );

  const auto& searchRange = CyclicPolygons::SearchRange(edgeLengths);
  BOOST_CHECK(std::fabs(searchRange.guess - 1 / (26.385 * 26.385)) < 1e-2);
}

BOOST_AUTO_TEST_CASE(automaticFindingWorks) {
  std::vector<double> edgeLengths {29, 30, 31, 32, 33};

  BOOST_CHECK(
    CyclicPolygons::maximumPentagonCircumradius(edgeLengths)
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
