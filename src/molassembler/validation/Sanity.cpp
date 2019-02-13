/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_TEST_MODULE SanityTests
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "chemical_symmetries/Symmetries.h"

#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include "temple/Functional.h"
#include "temple/Random.h"

using namespace std::string_literals;
using namespace Scine;
using namespace molassembler;

// Shows atom stereopermutator differences only
void explainDifference(
  const StereopermutatorList& a,
  const StereopermutatorList& b
) {
  std::cout << "First:" << std::endl;
  for(const auto& stereopermutator : a.atomStereopermutators()) {
    std::cout << stereopermutator.info() << "\n";
  }

  std::cout << "Second:" << std::endl;
  for(const auto& stereopermutator : b.atomStereopermutators()) {
    std::cout << stereopermutator.info() << "\n";
  }
  std::cout << "\n";
}

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

/* Test whether generating coordinates from a simple molecule and then
 * recovering all the stereopermutator data from the positions alone yields the
 * same StereopermutatorList as you started out with
 */
BOOST_AUTO_TEST_CASE( createPositionsAndFitNewMoleculeEqual ) {
  DistanceGeometry::Configuration DGConfiguration;
  DGConfiguration.partiality = DistanceGeometry::Partiality::All;

  for(const auto& symmetryName: Symmetry::allNames) {
    // Build an abstract asymmetric molecule (all ligands different) for the current molecule
    Molecule molecule(
      Utils::ElementType::Ru,
      Utils::ElementType::H,
      BondType::Single
    );

    for(unsigned i = 0; molecule.graph().N() - 1 < Symmetry::size(symmetryName); ++i) {
      molecule.addAtom(
        elements.at(i),
        0,
        BondType::Single
      );
    }

    if(!molecule.stereopermutators().empty()) {
      auto centralStereopermutator = molecule.stereopermutators().option(0);

      // Create a full list of the possible assignments
      std::vector<unsigned> assignments;
      assignments.resize(centralStereopermutator -> numAssignments());
      std::iota(
        assignments.begin(),
        assignments.end(),
        0
      );

      // Randomize
      temple::random::shuffle(assignments, randomnessEngine());

      /* Limit the number of assignments we're testing per symmetry to 10.
       * Otherwise, with maximally asymmetric square antiprismatic (5040),
       * we're never going to get done.
       *
       * Also because fitting stereopermutators to positions becomes ridiculously
       * expensive since it brute forces all assignments...
       *
       * Let's break it down to the number of Symmetry angle function calls
       * 10 assignments
       * 1000 Molecules
       * 1000 fit calls
       * 1000 * 5040 assignments tested
       * 1000 * 5040 * 8 * 7 ~= 3e8 angle function calls
       */
      if(centralStereopermutator -> numAssignments() > 100) {
        assignments.resize(10);
      }

      for(const auto& assignment : assignments) {
        molecule.assignStereopermutator(0, assignment);

        // For each possible arrangement of these ligands
        /* Create an ensemble of 3D positions using DG
         * and uniform distance setting
         */
        auto ensemble = DistanceGeometry::run(
          molecule,
          100,
          DGConfiguration
        );

        /* Check that for every PositionCollection, inferring the StereopermutatorList
         * from the generated coordinates yields the same StereopermutatorList you
         * started out with.
         */
        BOOST_CHECK_MESSAGE(
          temple::all_of(
            ensemble,
            [&](const auto& positionResult) -> bool {
              if(!positionResult) {
                std::cout << "Failed to generate a conformer: " << positionResult.error().message() << "\n";
              }

              auto inferredStereopermutatorList = molecule.inferStereopermutatorsFromPositions(positionResult.value());

              bool pass = molecule.stereopermutators() == inferredStereopermutatorList;

              if(!pass) {
                explainDifference(
                  molecule.stereopermutators(),
                  inferredStereopermutatorList
                );
              }

              return pass;
            }
          ),
          "Some reinterpretations of generated conformers did not yield identical stereopermutations!"
        );
      }
    }
  }
}
