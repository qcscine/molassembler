#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AdjacencyListTests
#include <boost/test/unit_test.hpp>

#include "AdjacencyList.h"
#include "template_magic/templateMagic.h"
#include "RNG/RNG.h"

#include <iostream>
#include <fstream>
#include <string>

/* To convert everything to RC checks, update boost
 * #include "rapidcheck/include/rapidcheck.h"
 * #include "rapidcheck/extras/boost_test/include/rapidcheck/boost_test.h"
 */

/* TODO
 * - Convert fixture tests to rapidcheck fixture properties so they can be 
 *   repeatedly checked randomly
 * - Consider separation of validation checking from the implemented class
 *   Should validate() and its properties be specified outside the class or be 
 *   expressed as pre- and postconditions for function calls?
 */

using namespace MoleculeManip;
using namespace std::string_literals;

namespace MoleculeManip {

struct AdjacencyListValidator {
  static boost::optional<std::string> everyAdjacencyHasReverse(
    const AdjacencyList& a
  ) {
    for(unsigned i = 0; i < a.size(); i++) {
      auto adjacents = a.getAdjacencies(i);
      for(const auto& adjacency: adjacents) {
        if(adjacency >= a.size()) {
          return "Listed adjacency is out of bounds"s;
        }
        if(!a.isAdjacent(adjacency, i)) {
          return "No inverse for the listed adjacency"s;
        }
      }
    }

    return boost::none;
  }

  // Empty rows can occur due to absolutely zero fault of the AdjacencyList if a
  // disconnect occurs, it is not inherently an invalid state.
  static boost::optional<std::string> noEmptyRows(
    const AdjacencyList& a
  ) {
    if(
      TemplateMagic::any_of(
        TemplateMagic::map(
          a._adjacencies,
          [](const std::vector<AtomIndexType>& localList) -> bool {
            return localList.size() == 0;
          }
        )
      )
    ) {
      return "There are empty rows in the ragged array!"s;
    } else {
      return boost::none;
    }
  }

  static std::vector<
    std::function<
      boost::optional<std::string>(
        const AdjacencyList&
      )
    >
  > validators; 

  static bool validate(const AdjacencyList& a) {
    for(const auto& validator : validators) {
      auto errorOption = validator(a);
      if(errorOption) {
        std::cout << errorOption.value() << std::endl;
        std::ofstream failureFile("validationFailure.dot");
        failureFile << a.dumpGraphviz() << std::endl;
        failureFile.close();
        return false;
      }
    }

    return true;
  }
};

std::vector<
  std::function<
    boost::optional<std::string>(
      const AdjacencyList&
    )
  >
> AdjacencyListValidator::validators {
  AdjacencyListValidator::everyAdjacencyHasReverse,
  AdjacencyListValidator::noEmptyRows
};

} // eo namespace MoleculeManip

// Create a random connected AdjacencyList
struct ALFixture {
  AdjacencyList adjacencies;

  ALFixture() {
    // generated adjacency lists parameters
    const AtomIndexType atomsLimit = 100;
    const unsigned cyclesLimit = 5;
    const unsigned edgesLimit = 6; // maximum edges per vertex number

    // start with two connected vertices
    AtomIndexType N = 2;
    adjacencies.resize(N);
    adjacencies.addAdjacency(0, 1);

    // extend at random with upper bound of 6 edges per vertex
    while(N < atomsLimit) {
      // select a random atom
      AtomIndexType selection = RNG::rng.getSingle<double>(0u, N - 1);

      // ensure less than 6 edges
      if(adjacencies.getAdjacencies(selection).size() >= edgesLimit) {
        // do not extend here
        continue;
      }

      // add a new slot
      auto newIndex = adjacencies.addSlot();
      
      // connect to selected random atom
      adjacencies.addAdjacency(selection, newIndex);

      N += 1;
    }

    // add some cycles
    unsigned nCycles = 0;
    while(nCycles < cyclesLimit) {
      // select two random atoms
      AtomIndexType i = RNG::rng.getSingle<double>(0u, N - 1);
      AtomIndexType j = RNG::rng.getSingle<double>(0u, N - 1);

      /* cannot connect
       * - equal indices 
       * - already adjacent ones 
       * - if either has maximum edge number
       */
      if(
        i == j
        || adjacencies.isAdjacent(i, j)
        || adjacencies.getAdjacencies(i).size() >= edgesLimit
        || adjacencies.getAdjacencies(j).size() >= edgesLimit
      ) {
        // skip a pair where any applies
        continue;
      }
      
      adjacencies.addAdjacency(i, j);
      nCycles += 1;
    }
  }
  ~ALFixture() {} // nothing to do
};

BOOST_FIXTURE_TEST_CASE(validFixture, ALFixture) {
  BOOST_CHECK(AdjacencyListValidator::validate(adjacencies));
}

// Non-trivial tests
BOOST_FIXTURE_TEST_CASE(indexInvalidation, ALFixture) {
  // generate a list of terminal vertices
  std::vector<AtomIndexType> terminalVertices;
  
  for(AtomIndexType i = 0; i < adjacencies.size(); i++) {
    if(adjacencies.getAdjacencies(i).size() == 1) {
      terminalVertices.push_back(i);
    }
  }

  if(terminalVertices.size() > 0) {
    // select a random terminal vertex
    auto selection = terminalVertices.at(
      static_cast<unsigned>(
        RNG::rng.getSingle<double>(0, terminalVertices.size() - 1)
      )
    );

    // remove all adjacencies involving it
    auto bondedIndices = adjacencies.getAdjacencies(selection);
    for(const auto& bondedIndex : bondedIndices) {
      adjacencies.removeAdjacency(selection, bondedIndex);
    }

    // minimalize
    adjacencies.indexInvalidationUpdate(selection);

    // test validity
    BOOST_CHECK(AdjacencyListValidator::validate(adjacencies));
  } 
}
