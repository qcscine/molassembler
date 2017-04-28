#define BOOST_TEST_MODULE AdjacencyListTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "AdjacencyList.h"
#include "template_magic/TemplateMagic.h"
#include "template_magic/Random.h"

// Temp, testing
#include <boost/graph/breadth_first_search.hpp>
#include "RepeatedElementCollection.h"
#include "Tree.h"
#include <algorithm>

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
    for(unsigned i = 0; i < a.numAtoms(); i++) {
      auto adjacents = a.getAdjacencies(i);
      for(const auto& adjacency: adjacents) {
        if(adjacency >= a.numAtoms()) {
          return "Listed adjacency is out of bounds"s;
        }
        if(!a.isAdjacent(adjacency, i)) {
          return "No inverse for the listed adjacency"s;
        }
      }
    }

    return boost::none;
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
        a.dumpGraphviz("validationFailure.dot");
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
};

} // namespace MoleculeManip

// Create a random connected AdjacencyList
struct ALFixture {
  AdjacencyList adjacencies;

  ALFixture() {
    // generated adjacency lists parameters
    const AtomIndexType atomsLimit = 10; // was 100
    const unsigned cyclesLimit = 1; // was 5
    const unsigned edgesLimit = 3; // was 6, maximum edges per vertex number

    // start with two connected vertices
    AtomIndexType N = 2;
    adjacencies.addAtom(Delib::ElementType::H);
    adjacencies.addAtom(Delib::ElementType::H);
    adjacencies.addBond(0, 1, BondType::Single);

    // extend at random with upper bound of 6 edges per vertex
    while(N < atomsLimit) {
      // select a random atom
      AtomIndexType selection = TemplateMagic::random.getSingle<double>(0u, N - 1);

      // ensure less than 6 edges
      if(adjacencies.getAdjacencies(selection).size() >= edgesLimit) {
        // do not extend here
        continue;
      }

      // add a new slot
      auto newIndex = adjacencies.addAtom(Delib::ElementType::H);
      
      // connect to selected random atom
      adjacencies.addBond(selection, newIndex, BondType::Single);

      N += 1;
    }

    // add some cycles
    unsigned nCycles = 0;
    while(nCycles < cyclesLimit) {
      // select two random atoms
      AtomIndexType i = TemplateMagic::random.getSingle<double>(0u, N - 1);
      AtomIndexType j = TemplateMagic::random.getSingle<double>(0u, N - 1);

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
      
      adjacencies.addBond(i, j, BondType::Single);
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
  
  for(AtomIndexType i = 0; i < adjacencies.numAtoms(); i++) {
    if(adjacencies.getAdjacencies(i).size() == 1) {
      terminalVertices.push_back(i);
    }
  }

  if(!terminalVertices.empty()) {
    // select a random terminal vertex
    auto selection = terminalVertices.at(
      static_cast<unsigned>(
        TemplateMagic::random.getSingle<double>(0, terminalVertices.size() - 1)
      )
    );

    // remove all adjacencies involving it
    auto bondedIndices = adjacencies.getAdjacencies(selection);
    for(const auto& bondedIndex : bondedIndices) {
      adjacencies.removeBond(selection, bondedIndex);
    }

    // test validity
    BOOST_CHECK(AdjacencyListValidator::validate(adjacencies));
  } 
}
