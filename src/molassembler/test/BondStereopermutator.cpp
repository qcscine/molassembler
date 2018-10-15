// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#define BOOST_TEST_MODULE BondStereopermutatorTestModule
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"

#include "molassembler/Conformers.h"
#include "molassembler/IO.h"
#include "molassembler/Molecule.h"
#include "molassembler/StereopermutatorList.h"

#include <iostream>

using BondIndex = molassembler::BondIndex;

struct Expectation {
  BondIndex edge;
  unsigned numPermutations;
  unsigned fittedAssignment;

  Expectation(
    BondIndex passEdge,
    unsigned permutations,
    unsigned assignment
  ) : edge(passEdge),
      numPermutations(permutations),
      fittedAssignment(assignment)
  {
    assert(permutations > 1);
  }

  Expectation(
    BondIndex passEdge,
    unsigned permutations
  ) : edge(passEdge),
      numPermutations(permutations)
  {
    assert(permutations == 1);
    // Define the value, but no comparisons with it should be performed
    fittedAssignment = std::numeric_limits<unsigned>::max();
  }
};

/* This is the current interpretation of yielded indices of permutations of
 * BondStereopermutator for all combinations of the symmetries trigonal planar and
 * bent.
 */
constexpr unsigned Z = 1;
constexpr unsigned E = 0;
constexpr unsigned stereogenic = 2;
constexpr unsigned nonStereogenic = 1;

std::map<std::string, Expectation> recognitionExpectations {
  {
    "but-2E-ene",
    {{0, 1}, stereogenic, E}
  },
  {
    "but-2Z-ene",
    {{0, 1}, stereogenic, Z}
  },
  {
    "E-diazene",
    {{0, 1}, stereogenic, E}
  },
  {
    "ethanimine",
    {{0, 2}, stereogenic, E}
  },
  {
    "ethene",
    {{0, 1}, nonStereogenic}
  },
  {
    "methanimine",
    {{0, 1}, nonStereogenic}
  }
  // formaldehyde is omitted, since there cannot be a BondStereopermutator on it
};

void checkExpectations(const boost::filesystem::path& filePath) {
  using namespace molassembler;
  using namespace std::string_literals;

  std::string moleculeName = filePath.stem().string();

  std::cout << "Processing " << moleculeName << std::endl;

  // Read the file
  auto mol = IO::read(filePath.string());

  // Check if expectations are met
  auto findIter = recognitionExpectations.find(filePath.stem().string());
  if(findIter == std::end(recognitionExpectations)) {
    // There should be no BondStereopermutator in this molecule
    auto bondStereopermutatorRange = mol.stereopermutators().bondStereopermutators();
    BOOST_CHECK(
      std::distance(
        std::begin(bondStereopermutatorRange),
        std::end(bondStereopermutatorRange)
      ) == 0
    );

    /* No DG work is needed since it doesn't involve BondStereopermutator (there
     * are none in the molecule)
     */
  } else {
    const Expectation& expectation = findIter->second;
    auto bondStereopermutatorOption = mol.stereopermutators().option(expectation.edge);

    BOOST_REQUIRE_MESSAGE(
      bondStereopermutatorOption,
      "There is no BondStereopermutator on the expected edge for " << moleculeName
    );

    BOOST_REQUIRE_MESSAGE(
      bondStereopermutatorOption->numStereopermutations() == expectation.numPermutations,
      "The expected number of permutations was not met for " << moleculeName
        << ": expected " << expectation.numPermutations << ", got "
        << bondStereopermutatorOption->numStereopermutations()
    );

    if(expectation.numPermutations == stereogenic) {
      auto assignmentOptional = bondStereopermutatorOption->assigned();

      BOOST_REQUIRE(assignmentOptional);

      BOOST_REQUIRE_MESSAGE(
        assignmentOptional.value() == expectation.fittedAssignment,
        "The stereopermutator is not assigned the expected value for "
          << moleculeName << ": Expected " << expectation.fittedAssignment
          << ", got " << *assignmentOptional << " instead"
      );
    }

    // Generate a conformation
    auto positionsResult = generateConformation(mol);

    // If DG fails, we're screwed
    if(!positionsResult) {
      BOOST_FAIL(positionsResult.error().message());
    }

    auto positions = positionsResult.value();

    // Reinterpret the molecule from the existing graph and the generated positions
    Molecule reinterpreted {
      mol.graph(),
      AngstromWrapper {positionsResult.value()}
    };

    bool pass = reinterpreted == mol;
    if(!pass) {
      std::cout << "Initial molecule: " << mol << "\nReinterpreted:"
        << reinterpreted;
    }
    BOOST_CHECK_MESSAGE(
      pass,
      "The reinterpreted molecule does not equal the initial molecule for "
        << moleculeName
    );
  }
}

BOOST_AUTO_TEST_CASE(stereopermutatorExpectationTests) {
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("ez_stereocenters")
  ) {
    BOOST_CHECK_NO_THROW(
      checkExpectations(currentFilePath)
    );
  }
}
