#include "CyclicPolygons.h"

#include <fstream>
#include <iostream>
#include "template_magic/Random.h"

int main() {

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
      const double geometricLimit = TemplateMagic::sum(edgeLengths);

      edgeLengths.emplace_back(
        TemplateMagic::random.getSingle<double>(
          lowerLimit,
          std::min(upperLimit, geometricLimit)
        )
      );
    }
  
    CyclicPolygons::analysis::writeSvrtanAnalysisFiles(
      edgeLengths,
      "scan-svrtan-"s + std::to_string(nTest)
    );
    CyclicPolygons::analysis::writeAngleAnalysisFiles(
      edgeLengths,
      "scan-angles-"s + std::to_string(nTest)
    );
  }


  /* Determine failure ratio of pentagon circumradius as a function of the
   * variance of a truncated normal distribution.
   */

  std::ofstream svrtanVarianceFile("edge-lengths-variance-success.csv");

  const double mean = (upperLimit - lowerLimit) / 2;
  const double lowerStddev = 0.01 * mean;
  const double upperStddev = 4 * mean;
  const unsigned nStddevSteps = 80;
  const unsigned nSamples = 1000;

  const double stepLength = (upperStddev - lowerStddev) / nStddevSteps;
  for(unsigned stepNumber = 0; stepNumber <= nStddevSteps; stepNumber++) {
    double currentStddev = lowerStddev + stepNumber * stepLength;

    std::normal_distribution<double> normalDistribution(mean, currentStddev);

    auto sampleTruncatedNormal = [&]() -> double {
      double sample;
      while(true) {
        sample = normalDistribution(TemplateMagic::random.randomEngine);
        if(lowerLimit <= sample && sample <= upperLimit) {
          return sample;
        }
      }
    };

    double successes = 0;
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

      auto circumradiusOption = CyclicPolygons::Pentagon::convexCircumradiusSvrtan(
        edgeLengths
      );

      if(circumradiusOption) {
        successes += 1;
      }
    }

    const double stddevFraction = currentStddev / mean;
    const double percentageCorrect = 100 * static_cast<double>(successes) / nSamples;

    // Write to file
    svrtanVarianceFile << stddevFraction << ", " << percentageCorrect << std::endl;
  }

  svrtanVarianceFile.close();
}
