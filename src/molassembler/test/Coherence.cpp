/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include "molassembler/IO.h"

#include "Utils/Geometry/ElementTypes.h"
#include "shapes/Data.h"

#include "temple/Functional.h"
#include "temple/Random.h"

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

// NOTE: This tests reflects upon a lot of things.
BOOST_AUTO_TEST_CASE(AtomStereopermutationCGFitCoherence) {
  using namespace Scine;
  DistanceGeometry::Configuration DGConfiguration;
  DGConfiguration.partiality = DistanceGeometry::Partiality::All;

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

    for(unsigned i = 0; molecule.graph().N() - 1 < Shapes::size(shape); ++i) {
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
    if(centralStereopermutator -> numAssignments() > 10) {
      assignments.resize(10);
    }

    for(const auto& assignment : assignments) {
      molecule.assignStereopermutator(0, assignment);

      // For each possible arrangement of these ligands, create an ensemble
      auto ensemble = DistanceGeometry::run(
        molecule,
        ensembleSize,
        DGConfiguration
      );

      /* Check that for every PositionCollection, inferring the StereopermutatorList
       * from the generated coordinates yields the same StereopermutatorList you
       * started out with.
       */
      unsigned matches = 0;
      unsigned mismatches = 0;
      boost::optional<AngstromWrapper> matchingPositions;
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
          + "_match.mol"
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
}
