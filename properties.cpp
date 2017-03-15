#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

#include "GenerateUniques.h"
#include "LogicalOperatorTests.h"

using namespace UniqueAssignments;

std::vector<char> createLowestSymmetryCharSequence(const unsigned& length) {
  std::vector<char> sequence {'A'};

  while(sequence.size() < length) {
    sequence.push_back(
      sequence.back() + 1
    );
  }

  return sequence;
}

int main() {
  for(const auto& symmetryName : Symmetry::allNames) {
    Assignment lowestSymmetryAssignment(
      symmetryName,
      createLowestSymmetryCharSequence(
        Symmetry::size(symmetryName)
      ),
      {} // No connections
    );

    auto uniques = UniqueAssignments::uniqueAssignments(
      lowestSymmetryAssignment
    );

    std::cout << Symmetry::name(symmetryName) << ": " << uniques.size() << std::endl;
  }

  return 0;
}
