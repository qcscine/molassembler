#define BOOST_TEST_MODULE CyclicPolygonsTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Minimal.h"
#include <iostream>
#include <fstream>
#include "temple/VectorView.h"
#include "temple/Random.h"
#include "temple/Stringify.h"

using namespace std::string_literals;

void writeAngleAnalysisFiles(
  const std::vector<double>& edgeLengths,
  const std::string& baseName
) {
  using namespace CyclicPolygons;

  const double minR = temple::max(edgeLengths) / 2 + 1e-10;

  /*const double lowerBound = std::max(
    detail::regularCircumradius(
      edgeLengths.size(),
      temple::min(edgeLengths)
    ),
    minR
  );*/
  const double lowerBound = minR;

  const double upperBound = std::max(
    detail::regularCircumradius(
      edgeLengths.size(),
      temple::max(edgeLengths)
    ), 
    minR
  );

  double rootGuess = detail::regularCircumradius(
    edgeLengths.size(),
    std::max(
      temple::average(edgeLengths),
      minR
    )
  );

  if(rootGuess < lowerBound) {
    rootGuess = lowerBound;
  } else if(rootGuess > upperBound) {
    rootGuess = upperBound;
  }

  auto rootSearchLambda = [&](const double& circumradius) -> std::tuple<double, double, double> {
    return std::make_tuple<double, double, double>(
      detail::centralAnglesDeviation(circumradius, edgeLengths),
      detail::centralAnglesDeviationDerivative(circumradius, edgeLengths),
      detail::centralAnglesDeviationSecondDerivative(circumradius, edgeLengths)
    );
  };

  const unsigned maxIter = 1000;
  boost::uintmax_t iterCount = maxIter;

  auto root = boost::math::tools::schroder_iterate(
    rootSearchLambda,
    rootGuess,
    lowerBound,
    upperBound,
    32, // bits precision
    iterCount
  );

  if(iterCount == static_cast<boost::uintmax_t>(maxIter)) {
    throw std::logic_error("Could not find polygon circumradius!");
  }

  std::ofstream scanFile(baseName + ".csv"s);
  scanFile << std::fixed << std::setprecision(8);

  const unsigned nScanSteps = 1000;
  const double stepSize = (upperBound - lowerBound) / nScanSteps;
  for(unsigned i = 0; i <= nScanSteps; i++) {
    const double currentR = lowerBound + i * stepSize;
    scanFile << currentR << ", " << detail::centralAnglesDeviation(currentR, edgeLengths) << std::endl;

  }

  scanFile.close();

  std::ofstream metaFile(baseName + "-meta.csv"s);
  metaFile << std::fixed << std::setprecision(8);
  for(unsigned i = 0; i < edgeLengths.size(); i++) {
    metaFile << edgeLengths[i];
    if(i != edgeLengths.size() - 1) {
      metaFile << ", ";
    }
  }
  metaFile << std::endl;
  metaFile << rootGuess << ", " << root << std::endl;
  metaFile.close();
}

BOOST_AUTO_TEST_CASE(centralAngleRootFinding) {
  unsigned failureIndex = 0;

  const double upperLimit = 5.6; // Fr-Fr single
  const double lowerLimit = 0.7; // H-H single

  const unsigned nTests = 100;

  for(unsigned nSides = 5; nSides < 10; ++nSides) {
    for(unsigned testNumber = 0; testNumber < nTests; testNumber++) {
      std::vector<double> edgeLengths = temple::random.getN<double>(lowerLimit, upperLimit, nSides);
      while(!CyclicPolygons::exists(edgeLengths)) {
        edgeLengths = temple::random.getN<double>(lowerLimit, upperLimit, nSides);
      }

      double circumradius;
      auto assignCircumradius = [&]() -> void {
        circumradius = CyclicPolygons::detail::convexCircumradius(
          edgeLengths
        );
      };

      BOOST_CHECK_NO_THROW(assignCircumradius());
      BOOST_CHECK(!std::isnan(circumradius));

      auto centralAnglesDeviation = CyclicPolygons::detail::centralAnglesDeviation(circumradius, edgeLengths);

      bool pass = std::fabs(centralAnglesDeviation) < 1e-6;

      BOOST_CHECK_MESSAGE(
        pass,
        "Central angles deviation for edge lengths " << temple::stringify(edgeLengths)
          << " is " << centralAnglesDeviation << ", whose norm is not less than 1e-6."
      );

      auto internalAngleSumDeviation = temple::sum(
        CyclicPolygons::detail::generalizedInternalAngles(edgeLengths, circumradius)
      ) - (nSides - 2) * M_PI;

      pass = pass && std::fabs(internalAngleSumDeviation < 1e-6);

      BOOST_CHECK_MESSAGE(
        std::fabs(internalAngleSumDeviation) < 1e-6,
        "Internal angle sum deviation from " << (nSides - 2) 
          <<  "π for edge lengths " << temple::stringify(edgeLengths)
          << " is " << internalAngleSumDeviation 
          << ", whose norm is not less than 1e-6"
      );

      if(!pass) {
        writeAngleAnalysisFiles(
          edgeLengths,
          "angle-failure-"s + std::to_string(failureIndex)
        );

        ++failureIndex;
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(specificCases) {
  using namespace CyclicPolygons;

  std::vector<double> edgeLengths = {3.627, 1.4, 1.533, 1.533, 1.4};

  BOOST_CHECK(exists(edgeLengths));

  auto circumradius = CyclicPolygons::detail::convexCircumradius(edgeLengths);

  std::cout << "r: " << circumradius << std::endl;

  auto internalAngles = CyclicPolygons::detail::generalizedInternalAngles(edgeLengths, circumradius);

  std::cout << temple::stringify(internalAngles) << std::endl;
}
