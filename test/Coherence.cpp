/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Conformers.h"
#include "Molassembler/DistanceGeometry/ConformerGeneration.h"
#include "Molassembler/DistanceGeometry/Error.h"

#include "Molassembler/RankingInformation.h"
#include "Molassembler/IO.h"
#include "Molassembler/IO/SmilesParser.h"

#include "Utils/Geometry/ElementTypes.h"
#include "Molassembler/Shapes/Data.h"

#include "Molassembler/Temple/Optionals.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Random.h"

#include "Fixtures.h"
#include "Utils/Typenames.h"

#include "boost/optional/optional_io.hpp"
#include <iostream>

using namespace std::string_literals;
using namespace Scine;
using namespace Molassembler;

// NOTE: Shows atom stereopermutator differences only
void explainDifference(
  const StereopermutatorList& a,
  const StereopermutatorList& b
) {
  std::cout << "First: ";
  for(const auto& stereopermutator : a.atomStereopermutators()) {
    std::cout << stereopermutator.info() << "\n";
  }

  std::cout << "Second: ";
  for(const auto& stereopermutator : b.atomStereopermutators()) {
    std::cout << stereopermutator.info() << "\n";
  }
  std::cout << "\n";
}

/* NOTE: This tests reflects upon a lot of things. It's run in the low
 * temperature approximation so that trigonal bipyramid and pentagonal bipyramid
 * don't get an easy out since they have one assignment only.
 */
BOOST_FIXTURE_TEST_CASE(AtomStereopermutationCGFitCoherence, LowTemperatureFixture) {
  using namespace Scine;
  DistanceGeometry::Configuration DgConfiguration;
  DgConfiguration.partiality = DistanceGeometry::Partiality::All;

  const std::array<Utils::ElementType, 9> elements {
    Utils::ElementType::F,
    Utils::ElementType::Cl,
    Utils::ElementType::Br,
    Utils::ElementType::I,
    Utils::ElementType::N,
    Utils::ElementType::C,
    Utils::ElementType::O,
    Utils::ElementType::S,
    Utils::ElementType::P
  };

#ifdef NDEBUG
  constexpr unsigned shapeSizeLimit = 8;
#else
  constexpr unsigned shapeSizeLimit = 4;
#endif

  constexpr unsigned ensembleSize = 5;

  for(const auto& shape: Shapes::allShapes) {
    if(Shapes::size(shape) > shapeSizeLimit) {
      continue;
    }

    // Build an abstract asymmetric molecule (all ligands different) for the current molecule
    Molecule molecule(
      Utils::ElementType::Ru,
      Utils::ElementType::H,
      BondType::Single
    );

    for(unsigned i = 0; molecule.graph().V() - 1 < Shapes::size(shape); ++i) {
      molecule.addAtom(elements.at(i), 0, BondType::Single);
    }

    assert(!molecule.stereopermutators().empty());

    auto centralStereopermutator = molecule.stereopermutators().option(0);
    assert(centralStereopermutator);
    auto assignments = Temple::iota<unsigned>(centralStereopermutator->numAssignments());

    // Pick a random assignment
    const unsigned assignment = Temple::Random::pick(assignments, randomnessEngine());

    molecule.assignStereopermutator(0, assignment);

    // For each possible arrangement of these ligands, create an ensemble
    auto ensemble = DistanceGeometry::run(
      molecule,
      ensembleSize,
      DgConfiguration,
      boost::none
    );

    /* Check that for every PositionCollection, inferring the StereopermutatorList
     * from the generated coordinates yields the same StereopermutatorList you
     * started out with.
     */
    unsigned matches = 0;
    unsigned mismatches = 0;
    boost::optional<AngstromPositions> matchingPositions;
    for(const auto& positionResult : ensemble) {
      if(!positionResult) {
        std::cout << "Failed to generate a conformer: " << positionResult.error().message() << "\n";
        continue;
      }

      auto inferredStereopermutatorList = molecule.inferStereopermutatorsFromPositions(positionResult.value());

      bool pass = molecule.stereopermutators() == inferredStereopermutatorList;

      if(pass) {
        if(matchingPositions == boost::none) {
          matchingPositions = positionResult.value();
        }
        ++matches;
      } else {
        explainDifference(
          molecule.stereopermutators(),
          inferredStereopermutatorList
        );

        std::string filename = (
          "shape_" + std::to_string(Shapes::nameIndex(shape))
          + "_" + std::to_string(assignment)
          + "_mismatch_" + std::to_string(mismatches)
          + ".mol"
        );
        IO::write(filename, molecule, positionResult.value());
        ++mismatches;
      }
    }

    if(matches != ensembleSize && matchingPositions) {
      std::string filename = (
        "shape_" + std::to_string(Shapes::nameIndex(shape))
        + "_" + std::to_string(assignment)
        + "match_.mol"
      );
      IO::write(
        filename,
        molecule,
        matchingPositions.value()
      );
    }

    BOOST_CHECK_MESSAGE(
      matches == ensembleSize,
      "Expected " << ensembleSize << " matches for assignment " << assignment << " of shape " << name(shape) << ", got " << matches << " matches."
    );
  }
}

void testIdenticalReinterpret(const Molecule& mol, const AtomIndex checkPosition) {
  auto fetchAssignment = [&](const Molecule& molecule) {
    return Temple::Optionals::flatMap(
      molecule.stereopermutators().option(checkPosition),
      [](const auto& f) { return f.assigned(); }
    );
  };

  const auto expectedAssignment = fetchAssignment(mol);
  const unsigned attempts = 3;
  Result<Utils::PositionCollection> conf {DgError::UnknownException};
  for(unsigned attempt = 0; attempt < attempts; ++attempt) {
    conf = generateRandomConformation(mol);
    if(conf) {
      break;
    }
  }
  if(!conf) {
    BOOST_FAIL(conf.error());
  }

  const Molecule reinterpreted {
    mol.graph(),
    AngstromPositions {conf.value()}
  };
  const auto reinterpretedAssignment = fetchAssignment(reinterpreted);
  BOOST_CHECK_EQUAL(reinterpretedAssignment, expectedAssignment);
}

#ifdef NDEBUG
BOOST_AUTO_TEST_CASE(BidentateAssignmentRecognized, *boost::unit_test::label("Molassembler")) {
  const std::string pincer_smiles = "[Ir]12([H])(Cl)P(C(C)(C)(C))(C(C)(C)(C))CC(=CC=C3)C1=C3CP2(C(C)(C)(C))C(C)(C)C";
  // NOTE: set shape at 0 to trigonal bipyramid
  auto pincer = IO::Experimental::parseSmilesSingleMolecule(pincer_smiles);
  pincer.setShapeAtAtom(0, Shapes::Shape::TrigonalBipyramid);

  // Find an assignment in which the phoshorus atoms are trans-arranged
  boost::optional<unsigned> expectedAssignmentOption;
  BOOST_REQUIRE(pincer.stereopermutators().option(0));
  const auto& permutator = pincer.stereopermutators().option(0).value();

  // Which sites at the permutator are the phosphorus atoms?
  std::vector<SiteIndex> phosphorusSites;
  const auto& ranking = permutator.getRanking();
  for(SiteIndex i {0}; i < ranking.sites.size(); ++i) {
    if(
      ranking.sites.at(i).size() == 1
      && pincer.graph().elementType(ranking.sites.at(i).front()) == Utils::ElementType::P
    ) {
      phosphorusSites.push_back(i);
    }
  }
  BOOST_REQUIRE(phosphorusSites.size() == 2);

  // Find a trans-arranged assignment
  const unsigned A = permutator.numAssignments();
  for(unsigned i = 0; i < A; ++i) {
    pincer.assignStereopermutator(0, i);
    if(std::fabs(permutator.angle(phosphorusSites.front(), phosphorusSites.back()) - M_PI) < 1e-10) {
      expectedAssignmentOption = i;
      break;
    }
  }
  BOOST_REQUIRE(expectedAssignmentOption);

  testIdenticalReinterpret(pincer, 0);
}
#endif

BOOST_AUTO_TEST_CASE(Shipscrews, *boost::unit_test::label("Molassembler")) {
  const std::string shipscrew_smiles = "[Fe@OH1+3]123(OC(=O)C(=O)O1)(OC(=O)C(=O)O2)OC(=O)C(=O)O3";
  auto shipscrew = IO::Experimental::parseSmilesSingleMolecule(shipscrew_smiles);

  const auto numAssignments = Temple::Optionals::flatMap(
    shipscrew.stereopermutators().option(0),
    [](const auto& p) { return p.numAssignments(); }
  );
  const auto expectedNumAssignments = boost::optional<unsigned>(2);
  BOOST_REQUIRE(numAssignments == expectedNumAssignments);

  for(unsigned i = 0; i < 2; ++i) {
    shipscrew.assignStereopermutator(0, i);
    testIdenticalReinterpret(shipscrew, 0);
  }
}
