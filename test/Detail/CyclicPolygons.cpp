/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Detail/CyclicPolygons.h"

#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/constexpr/Jsf.h"

#include <iostream>
#include <fstream>

using namespace Scine::Molassembler;
using namespace std::string_literals;

extern Temple::Generator<> generator;

void writeAngleAnalysisFiles(
  const std::vector<double>& edgeLengths,
  const std::string& baseName
) {
  using namespace CyclicPolygons;

  const double longestEdge = Temple::max(edgeLengths);
  const double minR = longestEdge / 2 + 1e-10;
  const double lowerBound = minR;
  const double upperBound = std::max(
    Detail::regularCircumradius(edgeLengths.size(), longestEdge),
    minR
  );

  double rootGuess = Detail::regularCircumradius(
    edgeLengths.size(),
    std::max(Temple::average(edgeLengths), minR)
  );

  if(rootGuess < lowerBound) {
    rootGuess = lowerBound;
  } else if(rootGuess > upperBound) {
    rootGuess = upperBound;
  }

  double circumradius;
  bool circumcenterInside;

  std::tie(circumradius, circumcenterInside) = CyclicPolygons::Detail::convexCircumradius(edgeLengths);

  std::ofstream scanFile(baseName + ".csv"s);
  scanFile << std::fixed << std::setprecision(8);

  const unsigned nScanSteps = 1000;
  const double stepSize = (upperBound - lowerBound) / nScanSteps;
  for(unsigned i = 0; i <= nScanSteps; i++) {
    const double currentR = lowerBound + i * stepSize;
    scanFile << currentR << ", "
      << Detail::circumcenterInside::centralAnglesDeviation(currentR, edgeLengths)
      << ", "
      << Detail::circumcenterOutside::centralAnglesDeviation(currentR, edgeLengths, longestEdge)
      << std::endl;
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
  metaFile << rootGuess << ", " << circumradius << std::endl;
  metaFile.close();
}

BOOST_AUTO_TEST_CASE(centralAngleRootFinding) {
  unsigned failureIndex = 0;

  const double upperLimit = 5.6; // Fr-Fr single
  const double lowerLimit = 0.7; // H-H single

  const unsigned nTests = 1000;

  for(unsigned nSides = 3; nSides < 10; ++nSides) {
    for(unsigned testNumber = 0; testNumber < nTests; testNumber++) {
      std::vector<double> edgeLengths = Temple::Random::getN<double>(lowerLimit, upperLimit, nSides, generator.engine);
      while(!CyclicPolygons::exists(edgeLengths)) {
        edgeLengths = Temple::Random::getN<double>(lowerLimit, upperLimit, nSides, generator.engine);
      }

      double circumradius = std::nan("");
      bool circumcenterInside = false;

      auto assignCircumradius = [&]() -> void {
        std::tie(circumradius, circumcenterInside) = CyclicPolygons::Detail::convexCircumradius(
          edgeLengths
        );
      };

      BOOST_CHECK_NO_THROW(assignCircumradius());
      BOOST_CHECK(!std::isnan(circumradius));

      double deviation;
      if(circumcenterInside) {
        deviation = CyclicPolygons::Detail::circumcenterInside::centralAnglesDeviation(
          circumradius,
          edgeLengths
        );
      } else {
        deviation = CyclicPolygons::Detail::circumcenterOutside::centralAnglesDeviation(
          circumradius,
          edgeLengths,
          Temple::max(edgeLengths)
        );
      }

      bool pass = std::fabs(deviation) < 1e-5;

      BOOST_CHECK_MESSAGE(
        pass,
        "Central angle deviation norm is not smaller than 1e-5 for " << Temple::stringify(edgeLengths)
          << ", circumcenter is inside: " << circumcenterInside << ", deviation: " << deviation
      );

      auto internalAngleSumDeviation = Temple::sum(
        CyclicPolygons::Detail::generalizedInternalAngles(edgeLengths, circumradius, circumcenterInside)
      ) - (nSides - 2) * M_PI;

      pass = pass && (std::fabs(internalAngleSumDeviation) < 1e-5);

      BOOST_CHECK_MESSAGE(
        std::fabs(internalAngleSumDeviation) < 1e-5,
        "Internal angle sum deviation from " << (nSides - 2)
          <<  "π for edge lengths " << Temple::stringify(edgeLengths)
          << " is " << internalAngleSumDeviation
          << ", whose norm is not less than 1e-5"
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
