#define BOOST_TEST_MODULE SanityTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "DistanceGeometry/generateConformation.h"
#include "BoundsFromSymmetry.h"
#include "IterateStereocenterPermutations.h"

#include "template_magic/Containers.h"
#include "template_magic/Random.h"

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

void explainDifference(
  const StereocenterList& a,
  const StereocenterList& b
) {
  std::cout << "A:" << std::endl;
  for(const auto& stereocenterPtr : a) {
    std::cout << stereocenterPtr << std::endl;
  }

  std::cout << "B:" << std::endl;
  for(const auto& stereocenterPtr : b) {
    std::cout << stereocenterPtr << std::endl;
  }
  std::cout << std::endl;
}

/* Test whether generating coordinates from a simple molecule and then
 * recovering all the stereocenter data from the positions alone yields the
 * same StereocenterList as you started out with
 */
BOOST_AUTO_TEST_CASE( createPositionsAndFitNewMoleculeEqual ) {
  for(const auto& symmetryName: Symmetry::allNames) {
    // Get an asymmetric molecule (all ligands different) for the current molecule
    auto molecule = DGDBM::asymmetricMolecule(symmetryName);

    if(molecule.getStereocenterList().size() > 0) {
      auto centralStereocenter = molecule.getStereocenterList().at(0);

      // Create a full list of the possible assignments
      std::vector<unsigned> assignments;
      assignments.resize(centralStereocenter -> numAssignments());
      std::iota(
        assignments.begin(),
        assignments.end(),
        0
      );

      // Randomize
      std::shuffle(
        assignments.begin(),
        assignments.end(),
        TemplateMagic::random.randomEngine
      );

      /* Limit the number of assignments we're testing per symmetry to 10.
       * Otherwise, with maximally asymmetric square antiprismatic (5040),
       * we're never going to get done.
       *
       * Also because fitting stereocenters to positions becomes ridiculously
       * expensive since it brute forces all assignments...
       *
       * Let's break it down to the number of Symmetry angle function calls
       * 10 assignments
       * 1000 Molecules
       * 1000 fit calls
       * 1000 * 5040 assignments tested
       * 1000 * 5040 * 8 * 7 ~= 3e8 angle function calls
       */
      if(centralStereocenter -> numAssignments() > 100) {
        assignments.resize(10);
      } 

      for(const auto& assignment : assignments) {
        molecule.assignStereocenterAtAtom(0, assignment);

        // For each possible arrangement of these ligands
        /* Create an ensemble of 3D positions using DG
         * and uniform distance setting
         */
        auto ensembleResult = detail::runDistanceGeometry(
          molecule,
          100,
          Partiality::All,
          false, // no y-inversion trick
          MoleculeSpatialModel::DistanceMethod::Uniform
        );

        BOOST_REQUIRE_MESSAGE(ensembleResult, ensembleResult.error().message());

        /* Check that for every PositionCollection, inferring the StereocenterList
         * from the generated coordinates yields the same StereocenterList you 
         * started out with.
         */
        auto mapped = TemplateMagic::map(
          ensembleResult.value(),
          [&](const auto& positions) -> bool {
            auto inferredStereocenterList = molecule.inferStereocentersFromPositions(
              positions
            );

            bool pass = molecule.getStereocenterList() == inferredStereocenterList;

            if(!pass) {
              explainDifference(
                molecule.getStereocenterList(),
                inferredStereocenterList
              );
            }

            return pass;
          }
        );

        /* The test passes only if this is true for all PositionCollections
         * yielded by DG
         */
        bool testPass = TemplateMagic::all_of(mapped);

        if(!testPass) {
          auto pass = TemplateMagic::count(mapped, true);

          std::cout << "Test fails!" << std::endl
            << std::setw(8) << " " << " " << Symmetry::name(symmetryName) 
            << std::endl
            << std::setw(8) << std::to_string(pass)+ "/100"
            << " comparisons with inferred StereocenterList pass" << std::endl;

          std::cout << "StereocenterList has entries:" << std::endl;
          for(const auto& stereocenter: molecule.getStereocenterList()) {
            std::cout << stereocenter << std::endl;
          }
        }

        BOOST_CHECK(testPass);
      }
    }
  }
}
