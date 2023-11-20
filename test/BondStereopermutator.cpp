/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"

#include "Molassembler/Conformers.h"
#include "Molassembler/IO.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Graph.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/StereopermutatorList.h"

#include <iostream>

using namespace Scine;
using BondIndex = Molassembler::BondIndex;

struct Expectation {
  BondIndex edge;
  unsigned numAssignments;
  unsigned fittedAssignment;

  Expectation(
    BondIndex passEdge,
    unsigned assignments,
    unsigned assignment
  ) : edge(passEdge),
      numAssignments(assignments),
      fittedAssignment(assignment)
  {
    assert(assignments > 1);
  }

  Expectation(
    BondIndex passEdge,
    unsigned assignments
  ) : edge(passEdge),
      numAssignments(assignments)
  {
    assert(assignments == 1);
    // Define the value, but no comparisons with it should be performed
    fittedAssignment = std::numeric_limits<unsigned>::max();
  }
};

/* This is the current interpretation of yielded indices of permutations of
 * BondStereopermutator for all combinations of the shapes triangle and
 * bent.
 */
constexpr unsigned Z = 1;
constexpr unsigned E = 0;
constexpr unsigned stereogenic = 2;
constexpr unsigned nonStereogenic = 1;

const std::map<std::string, Expectation> recognitionExpectations {
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
  using namespace Molassembler;
  using namespace std::string_literals;

  std::string moleculeName = filePath.stem().string();

  // Read the file
  auto mol = IO::read(filePath.string());

  // Check if expectations are met
  auto findIter = recognitionExpectations.find(filePath.stem().string());
  if(findIter == std::end(recognitionExpectations)) {
    // There should be no BondStereopermutator in this molecule
    auto bondStereopermutatorRange = mol.stereopermutators().bondStereopermutators();
    BOOST_CHECK(
      std::distance(
        bondStereopermutatorRange.first,
        bondStereopermutatorRange.second
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
      bondStereopermutatorOption->numAssignments() == expectation.numAssignments,
      "The expected number of permutations was not met for " << moleculeName
        << ": expected " << expectation.numAssignments << ", got "
        << bondStereopermutatorOption->numAssignments()
    );

    if(expectation.numAssignments == stereogenic) {
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
    auto positionsResult = generateRandomConformation(mol);

    // If DG fails, we're screwed
    if(!positionsResult) {
      BOOST_FAIL(positionsResult.error().message());
    }

    // Reinterpret the molecule from the existing graph and the generated positions
    Molecule reinterpreted {
      mol.graph(),
      AngstromPositions {positionsResult.value()}
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

BOOST_AUTO_TEST_CASE(BondStereopermutatorConsistency, *boost::unit_test::label("Molassembler")) {
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("ez_stereocenters")
  ) {
    BOOST_CHECK_NO_THROW(
      checkExpectations(currentFilePath)
    );
  }
}

BOOST_AUTO_TEST_CASE(BondStatePropagation, *boost::unit_test::label("Molassembler")) {
  using namespace Molassembler;

  auto mol = IO::read("ez_stereocenters/but-2E-ene.mol");

  // Alter a hydrogen at the bond stereopermutator
  const StereopermutatorList& stereopermutators = mol.stereopermutators();

  auto bondStereopermutatorRange = stereopermutators.bondStereopermutators();
  BOOST_REQUIRE(std::distance(bondStereopermutatorRange.first, bondStereopermutatorRange.second) > 0);

  const BondStereopermutator& mainStereopermutator = *bondStereopermutatorRange.first;

  BOOST_REQUIRE(mainStereopermutator.assigned());

  unsigned priorAssignment = mainStereopermutator.assigned().value();

  // Pick a side
  const AtomIndex side = mainStereopermutator.placement().first;

  // Find a hydrogen substituent
  boost::optional<AtomIndex> hydrogenSubstituent;
  for(const AtomIndex substituent : mol.graph().adjacents(side)) {
    if(mol.graph().elementType(substituent) == Utils::ElementType::H) {
      hydrogenSubstituent = substituent;
      break;
    }
  }

  BOOST_REQUIRE(hydrogenSubstituent);

  // Replace the hydrogen substituent with a fluorine
  mol.setElementType(*hydrogenSubstituent, Utils::ElementType::F);

  // All references are, in principle, invalidated. Just being extra careful.
  auto postPermutatorRange = stereopermutators.bondStereopermutators();
  // The new stereopermutator must still be assigned, and have a different assignment
  BOOST_REQUIRE(std::distance(postPermutatorRange.first, postPermutatorRange.second) > 0);
  const BondStereopermutator& postPermutator = *postPermutatorRange.first;
  BOOST_REQUIRE(postPermutator.assigned());

  // In this particular case, we know that the final assignment has to be different
  BOOST_CHECK(postPermutator.assigned().value() != priorAssignment);
}

BOOST_AUTO_TEST_CASE(BondStatePropagationOnRemoval, *boost::unit_test::label("Molassembler")) {
  using namespace Molassembler;

  auto mol = IO::read("various/benzene.mol");
  mol.removeAtom(6);
  mol.removeAtom(0);
  mol.removeAtom(5);
  mol.removeAtom(0);
  mol.removeAtom(4);
  mol.removeAtom(0);
  mol.removeAtom(3);
  mol.removeAtom(0);
  mol.removeAtom(2);
  mol.removeAtom(0);
  mol.removeAtom(1);
}

BOOST_AUTO_TEST_CASE(BondStatePropagationOnRemoval2, *boost::unit_test::label("Molassembler")) {
  using namespace Molassembler;

  auto mol = IO::read("various/benzene.mol");
  mol.removeAtom(11);
  mol.removeAtom(10);
  mol.removeAtom(9);
  mol.removeAtom(8);
  mol.removeAtom(7);
  mol.removeAtom(6);
  mol.removeAtom(5);
  mol.removeAtom(4);
  mol.removeAtom(3);
  mol.removeAtom(2);
  mol.removeAtom(1);
}

BOOST_AUTO_TEST_CASE(StereocentersInSmallCycles, *boost::unit_test::label("Molassembler")) {
  // Flat map from cycle size to number of assignments
  const std::vector<unsigned> expectedAssignmentsMap {
    0, 0, 0, 1, 1, 1, 1, 2, 2
  };

  using namespace Molassembler;
  /* In small cycles, double bond stereopermutators should exclude the E
   * stereopermutation and have only a single assignment
   */
  for(unsigned cycleSize = 3; cycleSize <= 8; ++cycleSize) {
    Molecule mol {
      Utils::ElementType::C,
      Utils::ElementType::C,
      BondType::Double
    };

    // Add cycle atoms
    AtomIndex last = 0;
    for(unsigned i = 0; i < cycleSize - 2; ++i) {
      last = mol.addAtom(Utils::ElementType::C, last);
    }

    // Close the cycle
    mol.addBond(last, 1);

    // Set geometries
    for(unsigned i = 0; i < cycleSize; ++i) {
      mol.setShapeAtAtom(i, Shapes::Shape::Bent);
    }

    auto bondStereopermutatorOption = mol.stereopermutators().option(BondIndex {0, 1});

    BOOST_REQUIRE_MESSAGE(
      bondStereopermutatorOption,
      "Expected a BondStereopermutator on {0, 1}, got None"
    );

    BOOST_CHECK_MESSAGE(
      bondStereopermutatorOption->numAssignments() == expectedAssignmentsMap.at(cycleSize),
      "Expected " << expectedAssignmentsMap.at(cycleSize)
        << " assignments for cycle of size " << cycleSize
        << ", got " << bondStereopermutatorOption->numAssignments()
        << "instead."
    );
  }
}
