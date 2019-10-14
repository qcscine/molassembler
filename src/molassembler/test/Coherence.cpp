/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include "Utils/Geometry/ElementTypes.h"
#include "chemical_symmetries/Symmetries.h"

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

/* NOTE: This tests reflects upon a lot of things.
 */
BOOST_AUTO_TEST_CASE(CGReinterpretYieldsSameShapes) {
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

  for(const auto& shape: Symmetry::allShapes) {
    // Build an abstract asymmetric molecule (all ligands different) for the current molecule
    Molecule molecule(
      Utils::ElementType::Ru,
      Utils::ElementType::H,
      BondType::Single
    );

    for(unsigned i = 0; molecule.graph().N() - 1 < Symmetry::size(shape); ++i) {
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
        20,
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
