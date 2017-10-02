#include "Math.h"

#include "template_magic/Containers.h"
#include "BTree.h"

#include <fstream>
#include <iomanip>
#include <vector>

// TODO TEMP
#include <iostream>

std::vector<double> log10Spaced(
  const double& minimum,
  const double& maximum,
  const unsigned& nSteps
) {
  auto logMin = ConstexprMagic::Math::log10(minimum);
  auto logMax = ConstexprMagic::Math::log10(maximum);

  auto logStepLength = (logMax - logMin) / nSteps;

  std::vector<double> spaced;

  for(unsigned i = 0; i < nSteps; ++i) {
    spaced.push_back(
      std::pow(10, logMin + i * logStepLength)
    );
  }

  return spaced;
}

void writeAsinFile() {
  const double lowerTrigLimit = -1 + std::numeric_limits<double>::epsilon();
  const double upperTrigLimit = 1 - std::numeric_limits<double>::epsilon();
  const double stepSize = 1e-3;

  std::ofstream asinFile("asin.csv");

  asinFile << std::scientific << std::setprecision(3);

  for(double value = lowerTrigLimit; value <= upperTrigLimit; value += stepSize) {
    asinFile << value << "," 
      << ConstexprMagic::Math::asin(value) << ","
      << std::asin(value) << "\n";
  }

  asinFile.close();
}

void writeBTreeFile() {
  const size_t minMinDegree = 2;
  const size_t maxMinDegree = 20;
  const double minElements = 10;
  const double maxElements = 1e5;

  std::ofstream bTreeFile("btree.csv");

  bTreeFile << std::scientific << std::setprecision(14);

  for(size_t minDegree = minMinDegree; minDegree <= maxMinDegree; ++minDegree) {
    auto numElements = TemplateMagic::map(
      log10Spaced(
        minElements,
        maxElements,
        500
      ),
      [](const double& numElements) -> size_t {
        return ConstexprMagic::Math::floor(numElements);
      }
    );

    auto percentageUsed = TemplateMagic::map(
      numElements,
      [&](const size_t& requestedSize) -> double {
        auto spaceAllocated = ConstexprMagic::BTreeProperties::maxNodesInTree(
          ConstexprMagic::BTreeProperties::maxHeightBound(
            requestedSize,
            minDegree
          ),
          minDegree
        ) * (2 * minDegree - 1);

        return 100.0 * requestedSize / spaceAllocated;
      }
    );

    bTreeFile << TemplateMagic::condenseIterable(percentageUsed) << std::endl;
  }

  bTreeFile.close();
}

int main() {
  writeAsinFile();
  writeBTreeFile();
}
