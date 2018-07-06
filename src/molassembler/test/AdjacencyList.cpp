#define BOOST_TEST_MODULE MoleculeGraphTests
#include <boost/test/unit_test.hpp>

#include "Molecule.h"
#include "temple/Random.h"

// Temp, testing
#include <boost/graph/breadth_first_search.hpp>
#include "RepeatedElementCollection.h"
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

using namespace molassembler;
using namespace std::string_literals;

namespace molassembler {

struct MoleculeValidator {
  static boost::optional<std::string> everyAdjacencyHasReverse(
    const Molecule& a
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
        const Molecule&
      )
    >
  > validators;

  static bool validate(const Molecule& a) {
    for(const auto& validator : validators) {
      auto errorOption = validator(a);
      if(errorOption) {
        std::cout << errorOption.value() << std::endl;
        std::ofstream outFile ("validationFailure.dot");
        outFile << a.dumpGraphviz();
        outFile.close();
        return false;
      }
    }

    return true;
  }
};

std::vector<
  std::function<
    boost::optional<std::string>(
      const Molecule&
    )
  >
> MoleculeValidator::validators {
  MoleculeValidator::everyAdjacencyHasReverse,
};

} // namespace molassembler

// Create a random connected Molecule
struct ALFixture {
  Molecule molecule {
    Delib::ElementType::H,
    Delib::ElementType::H,
    BondType::Single
  };

  ALFixture() {
    // generated adjacency lists parameters
    const AtomIndexType atomsLimit = 10; // was 100
    const unsigned cyclesLimit = 1; // was 5
    const unsigned edgesLimit = 3; // was 6, maximum edges per vertex number

    // start with two connected vertices
    AtomIndexType N = 2;

    // extend at random with upper bound of 6 edges per vertex
    while(N < atomsLimit) {
      // select a random atom
      AtomIndexType selection = rng.getSingle<double>(0u, N - 1);

      // ensure less than 6 edges
      if(molecule.getAdjacencies(selection).size() >= edgesLimit) {
        // do not extend here
        continue;
      }

      // add a new atom connected to the random atom
      molecule.addAtom(
        Delib::ElementType::H,
        selection,
        BondType::Single
      );

      N += 1;
    }

    // add some cycles
    unsigned nCycles = 0;
    while(nCycles < cyclesLimit) {
      // select two random atoms
      AtomIndexType i = rng.getSingle<double>(0u, N - 1);
      AtomIndexType j = rng.getSingle<double>(0u, N - 1);

      /* cannot connect
       * - equal indices
       * - already adjacent ones
       * - if either has maximum edge number
       */
      if(
        i == j
        || molecule.isAdjacent(i, j)
        || molecule.getAdjacencies(i).size() >= edgesLimit
        || molecule.getAdjacencies(j).size() >= edgesLimit
      ) {
        // skip a pair where any applies
        continue;
      }

      molecule.addBond(i, j, BondType::Single);
      nCycles += 1;
    }
  }
  ~ALFixture() = default;
};

BOOST_FIXTURE_TEST_CASE(validFixture, ALFixture) {
  BOOST_CHECK(MoleculeValidator::validate(molecule));
}

// Non-trivial tests
BOOST_FIXTURE_TEST_CASE(indexInvalidation, ALFixture) {
  // generate a list of terminal vertices
  std::vector<AtomIndexType> terminalVertices;

  for(AtomIndexType i = 0; i < molecule.numAtoms(); i++) {
    if(molecule.getAdjacencies(i).size() == 1) {
      terminalVertices.push_back(i);
    }
  }

  if(!terminalVertices.empty()) {
    // select a random terminal vertex
    auto selection = terminalVertices.at(
      rng.getSingle<unsigned>(0, terminalVertices.size() - 1)
    );

    // remove all bonds involving it
    molecule.removeAtom(selection);

    // test validity
    BOOST_CHECK(MoleculeValidator::validate(molecule));
  }
}
