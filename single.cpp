#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

#include "GenerateUniques.h"
#include "SymmetryInformation.h"

#include "LogicalOperatorTests.h"

using namespace UniqueAssignments;

/* TODO
 * - investigate failures!
 */

template<typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& vector) {
  for(unsigned i = 0; i < vector.size(); i++) {
    os << vector[i];
    if(i != vector.size() - 1) {
      os << ", ";
    }
  }

  return os;
}

template<typename Iterable>
void printAll(const Iterable& container) {
  for(const auto& item : container) {
    std::cout << item << std::endl;
  }
}


template<class Symmetry>
void run_tests(
  const std::vector<
    std::tuple<
      std::vector<char>, // characters
      std::set<
        std::pair<unsigned, unsigned>
      >, // pairs
      unsigned // expectedUnique
    >
  >& test_cases
) {

  for(const auto& tuple: test_cases) {
    std::vector<char> characters;
    std::set<
      std::pair<unsigned, unsigned>
    > pairs;
    unsigned expectedUnique;

    std::tie(characters, pairs, expectedUnique) = tuple;

    // instantiate
    Assignment<Symmetry> assignment = (pairs.size() == 0)
      ? Assignment<Symmetry>(characters)
      : Assignment<Symmetry>(characters, pairs);

    auto unique = uniqueAssignments(assignment);
    BOOST_CHECK(unique.size() == expectedUnique );
    if(unique.size() != expectedUnique) {
      std::cout << "Mismatch: Expected " << expectedUnique
        << " assignments for: \n" << assignment << ", got " 
        << unique.size() << " assignments:" << std::endl;
      for(const auto& uniqueAssignment: unique) {
        std::cout << uniqueAssignment << std::endl;
      }
    } 
  }
}

BOOST_AUTO_TEST_CASE( columnSmallerConsistency ) {
  Assignment<PermSymmetry::Octahedral> single {
    {'A', 'A', 'A', 'A', 'A', 'A'},
    {
      std::make_pair(0,1),
      std::make_pair(2,3),
      std::make_pair(4,5)
    }
  };

  bool pass = true;

  do {
    for(unsigned i = 0; i < 6 && pass; i++) {
      for(unsigned j = i + 1; j < 6 && pass; j++) {
        if(!OperatorTests::XOR(
            single.columnSmaller(i, j),
            single.columnSmaller(j, i),
            !single.columnSmaller(i, j) && !single.columnSmaller(j, i)
        )) {
          pass = false;
        }
      }
    }
  } while(single.nextPermutation() && pass);

  BOOST_CHECK(pass);
}

BOOST_AUTO_TEST_CASE( lowestPermutation) {
  Assignment<PermSymmetry::Octahedral> single {
    {'A', 'A', 'A', 'A', 'A', 'A'},
    {
      std::make_pair(0,1),
      std::make_pair(2,3),
      std::make_pair(4,5)
    }
  };

  single.lowestPermutation();

  BOOST_CHECK(single.isSortedAsc());
}

BOOST_AUTO_TEST_CASE( octahedral_multidentate ) {
  run_tests<PermSymmetry::Octahedral>(
    {
      // M(A-B)_2 CD
      std::make_tuple(
        std::vector<char>({'A', 'B', 'A', 'B', 'C', 'D'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3)
        }),
        11 // TODO ERROR: get 10
      )
    }
  );
}
