#include "Math.h"

#include <fstream>
#include <iomanip>

int main() {
  const double lowerTrigLimit = -1 + std::numeric_limits<double>::epsilon();
  const double upperTrigLimit = 1 - std::numeric_limits<double>::epsilon();
  const double stepSize = 1e-3;

  std::ofstream asinFile("asin.csv");

  asinFile << std::scientific << std::setprecision(14);

  for(double value = lowerTrigLimit; value <= upperTrigLimit; value += stepSize) {
    asinFile << value << "," 
      << ConstexprMagic::Math::asin(value) << ","
      << std::asin(value) << "\n";
  }

  asinFile.close();
}
