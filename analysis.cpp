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
      const double geometricLimit = TemplateMagic::numeric::sum(edgeLengths);

      edgeLengths.emplace_back(
        TemplateMagic::random.getSingle<double>(
          lowerLimit,
          std::min(upperLimit, geometricLimit)
        )
      );
    }
  
    std::string baseString = "scan-random-"s + std::to_string(nTest);

    CyclicPolygons::analysis::writeAnalysisFiles(edgeLengths, baseString);
  }


  /* Determine failure ratio of pentagon circumradius as a function of the
   * variance of a truncated normal distribution.
   */
  
  // Randomness set-up
  std::vector<unsigned> seeds;
  std::mt19937 randomEngine;

  { // quick local scope to avoid namespace pollution
#ifdef NDEBUG
    std::random_device randomDevice;
    for(unsigned n = 0; n < 5; n++) seeds.emplace_back(randomDevice());
#else 
    seeds.emplace_back(2721813754);
#endif
    std::seed_seq seedSequence(seeds.begin(), seeds.end());
    randomEngine.seed(seedSequence);
  }

  std::ofstream varianceFile("edge-lengths-variance-success.csv");

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
        sample = normalDistribution(randomEngine);
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

      while(!CyclicPolygons::cyclicPolygonConstructible(edgeLengths)) {
        edgeLengths.clear();
        while(edgeLengths.size() < 5) {
          edgeLengths.emplace_back(sampleTruncatedNormal());
        }
      }

      auto circumradiusOption = CyclicPolygons::Pentagon::maximumCircumradius(
        edgeLengths
      );

      if(circumradiusOption) {
        successes += 1;
      }
    }

    const double stddevFraction = currentStddev / mean;
    const double percentageCorrect = 100 * static_cast<double>(successes) / nSamples;

    // Write to file
    varianceFile << stddevFraction << ", " << percentageCorrect << std::endl;
  }

  varianceFile.close();
}
