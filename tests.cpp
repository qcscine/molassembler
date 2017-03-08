#define BOOST_TEST_MODULE SymmetryTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Symmetries.h"
#include <set>
#include <iostream>
// for accumulate, include numeric
#include <numeric>

using namespace Symmetry;

std::vector<unsigned> rotate(
  const std::vector<unsigned>& toRotate,
  const std::vector<unsigned>& rotationVector
) {
  std::vector<unsigned> rotated (toRotate.size()); 

  for(unsigned i = 0; i < toRotate.size(); i++) {
    rotated[i] = toRotate[
      rotationVector[i]
    ];
  }

  return rotated;
}


BOOST_AUTO_TEST_CASE( symmetrySanityTests ) {
  // every rotation vector size must equal size of symmetry
  for(const auto& name : allNames) {
    for(const auto& rotationVector : rotations(name)) {
      BOOST_CHECK(rotationVector.size() == size(name));
    }
  }

  // every rotation may have every number 0 -> (size of symmetry - 1) only once
  for(const auto& name : allNames) {
    std::set<unsigned> members;
    for(unsigned i = 0; i < size(name); i++) {
      members.insert(members.end(), i);
    }

    for(const auto& rotationVector : rotations(name)) {
      std::set<unsigned> converted {
        rotationVector.begin(),
        rotationVector.end()
      };

      BOOST_CHECK(converted.size() == size(name)); // no duplicates
      
      BOOST_CHECK(
        std::accumulate(
          rotationVector.begin(),
          rotationVector.end(),
          true,
          [&members](const bool& carry, const unsigned& rotationElement) {
            return carry && members.count(rotationElement) == 1;
          }
        )
      );
    }
  }

  /* every rotation must return to the original after a finite number of
   * applications
   */
  unsigned maxIter = 100;
  for(const auto& name : allNames) {
    std::vector<unsigned> initialConfiguration (
      size(name),
      0
    );

    std::iota(
      initialConfiguration.begin(),
      initialConfiguration.end(),
      0
    );

    for(const auto& rotationVector : rotations(name)) {
      // copy in from initial
      auto configuration = initialConfiguration;

      bool pass = false;
      for(unsigned N = 0; N < maxIter; N++) {
        configuration = rotate(configuration, rotationVector);
        if(configuration == initialConfiguration) {
          pass = true;
          break;
        }
      }

      BOOST_CHECK(pass);
    }
  }

  // every angle function must be symmetrical on input of valid unsigned indices
  for(const auto& symmetryName: allNames) {
    bool passesAll = true;

    for(unsigned i = 0; i < size(symmetryName) && passesAll; i++) {
      for(unsigned j = i + 1; j < size(symmetryName); j++) {
        if(angleFunction(symmetryName)(i, j) != angleFunction(symmetryName)(j, i)) {
          passesAll = false;
          std::cout << name(symmetryName) 
            << " is not symmetrical w.r.t. input indices: falsified by (" 
            << i << ", " << j <<") -> (" << angleFunction(symmetryName)(i, j) 
            << ", " << angleFunction(symmetryName)(j, i) << ")." << std::endl;
          break;
        }
      }
    }

    BOOST_CHECK(passesAll);
  }

  // every angle function must return 0 for identical indices
  for(const auto& symmetryName: allNames) {
    bool passesAll = true;

    for(unsigned i = 0; i < size(symmetryName); i++) {
      if(angleFunction(symmetryName)(i, i) != 0) {
        passesAll = false;
        std::cout << name(symmetryName)
          << "'s angle function does not return zero for identical indices ("
          << i << ", " << i << ")." << std::endl;
        break;
      }
    }

    BOOST_CHECK(passesAll);
  }
}
