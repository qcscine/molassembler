/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Conformers.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include "molassembler/RankingInformation.h"
#include "molassembler/IO.h"
#include "molassembler/IO/SmilesParser.h"

#include "Utils/Geometry/ElementTypes.h"
#include "shapes/Data.h"

#include "temple/Optionals.h"
#include "temple/Functional.h"
#include "temple/Random.h"

#include "Fixtures.h"

#include "boost/optional/optional_io.hpp"
#include <iostream>

using namespace std::string_literals;
using namespace Scine;
using namespace molassembler;

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
  distance_geometry::Configuration DgConfiguration;
  DgConfiguration.partiality = distance_geometry::Partiality::All;

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

  constexpr unsigned ensembleSize = 10;

  for(const auto& shape: shapes::allShapes) {
    if(shapes::size(shape) > shapeSizeLimit) {
      continue;
    }

    // Build an abstract asymmetric molecule (all ligands different) for the current molecule
    Molecule molecule(
      Utils::ElementType::Ru,
      Utils::ElementType::H,
      BondType::Single
    );

    for(unsigned i = 0; molecule.graph().N() - 1 < shapes::size(shape); ++i) {
      molecule.addAtom(
        elements.at(i),
        0,
        BondType::Single
      );
    }

    assert(!molecule.stereopermutators().empty());

    auto centralStereopermutator = molecule.stereopermutators().option(0);
    assert(centralStereopermutator);
    auto assignments = temple::iota<unsigned>(centralStereopermutator->numAssignments());

    // Randomize
    temple::random::shuffle(assignments, randomnessEngine());

    /* Limit the number of assignments we're testing per shape to 10.
     * Otherwise, with maximally asymmetric square antiprismatic (5040),
     * we're never going to get done.
     *
     * Also because fitting stereopermutators to positions becomes ridiculously
     * expensive since it brute forces all assignments...
     */
    if(assignments.size() > 10) {
      assignments.resize(10);
    }

    for(const auto& assignment : assignments) {
      molecule.assignStereopermutator(0, assignment);

      // For each possible arrangement of these ligands, create an ensemble
      auto ensemble = distance_geometry::run(
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
            "shape_" + std::to_string(shapes::nameIndex(shape))
            + "_" + std::to_string(assignment)
            + "mismatch__" + std::to_string(mismatches)
            + ".mol"
          );
          io::write(filename, molecule, positionResult.value());
          ++mismatches;
        }
      }

      if(matches != ensembleSize && matchingPositions) {
        std::string filename = (
          "shape_" + std::to_string(shapes::nameIndex(shape))
          + "_" + std::to_string(assignment)
          + "match_.mol"
        );
        io::write(
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
}

void testIdenticalReinterpret(const Molecule& mol, const AtomIndex checkPosition) {
  auto fetchAssignment = [&](const Molecule& molecule) {
    return temple::optionals::flatMap(
      molecule.stereopermutators().option(checkPosition),
      [](const auto& f) { return f.assigned(); }
    );
  };

  const auto expectedAssignment = fetchAssignment(mol);
  if(auto conf = generateRandomConformation(mol)) {
    Molecule reinterpreted {
      mol.graph(),
      AngstromPositions {conf.value()}
    };
    auto reinterpretedAssignment = fetchAssignment(reinterpreted);
    BOOST_CHECK_EQUAL(reinterpretedAssignment, expectedAssignment);
  } else {
    BOOST_FAIL(conf.error());
  }
}

BOOST_AUTO_TEST_CASE(BidentateAssignmentRecognized) {
  const std::string pincer_smiles = "[Ir]12([H])(Cl)P(C(C)(C)(C))(C(C)(C)(C))CC(=CC=C3)C1=C3CP2(C(C)(C)(C))C(C)(C)C";
  // NOTE: set shape at 0 to trigonal bipyramid
  auto pincer = io::experimental::parseSmilesSingleMolecule(pincer_smiles);
  pincer.setShapeAtAtom(0, shapes::Shape::TrigonalBipyramid);

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

BOOST_AUTO_TEST_CASE(Shipscrews) {
  const std::string shipscrew_smiles = "[Fe@OH1+3]123(OC(=O)C(=O)O1)(OC(=O)C(=O)O2)OC(=O)C(=O)O3";
  auto shipscrew = io::experimental::parseSmilesSingleMolecule(shipscrew_smiles);

  const auto numAssignments = temple::optionals::flatMap(
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
