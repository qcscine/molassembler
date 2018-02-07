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

std::vector<char> twoIdenticalCharSequence(const unsigned& length) {
  std::vector<char> sequence {'A', 'A'};

  while(sequence.size() < length) {
    sequence.push_back(
      sequence.back() + 1
    );
  }

  return sequence;
}

int main() {
  std::cout << "All different characters. Number of unique assignments:" << std::endl;
  for(const auto& symmetryName : Symmetry::allNames) {
    Assignment lowestSymmetryAssignment(
      symmetryName,
      createLowestSymmetryCharSequence(
        Symmetry::size(symmetryName)
      ),
      {} // No connections
    );

    auto uniques = UniqueAssignments::uniqueAssignments(
      lowestSymmetryAssignment,
      symmetryName
    );

    std::cout << Symmetry::name(symmetryName) << ": " << uniques.size() << std::endl;
  }

  std::cout << std::endl << "Two chars identical, all others different. Number of unique assignments:" << std::endl;
  for(const auto& symmetryName : Symmetry::allNames) {
    Assignment lowestSymmetryAssignment(
      symmetryName,
      twoIdenticalCharSequence(
        Symmetry::size(symmetryName)
      ),
      {} // No connections
    );

    auto uniques = UniqueAssignments::uniqueAssignments(
      lowestSymmetryAssignment,
      symmetryName
    );

    std::cout << Symmetry::name(symmetryName) << ": " << uniques.size() << std::endl;
  }

  std::cout << std::endl << "Two linked chars identical, all others different. Number of unique assignments:" << std::endl;
  for(const auto& symmetryName : Symmetry::allNames) {
    Assignment lowestSymmetryAssignment(
      symmetryName,
      twoIdenticalCharSequence(
        Symmetry::size(symmetryName)
      ),
      {
        {0, 1}
      }
    );

    auto uniques = UniqueAssignments::uniqueAssignments(
      lowestSymmetryAssignment,
      symmetryName
    );

    std::cout << Symmetry::name(symmetryName) << ": " << uniques.size() << std::endl;
  }

  return 0;
}
