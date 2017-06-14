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
}
