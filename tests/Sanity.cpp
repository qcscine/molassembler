#define BOOST_TEST_MODULE DGRefinementProblemTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"

#include "BoundsFromSymmetry.h"

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

/* Test whether generating coordinates from a simple molecule and then
 * recovering all the stereocenter data from the positions alone yields the
 * same StereocenterList as you started out with
 */
BOOST_AUTO_TEST_CASE( createPositionsAndFitNewMoleculeEqual ) {
  for(const auto& symmetryName: Symmetry::allNames) {
    auto molecule = DGDBM::symmetricMolecule(symmetryName);

    auto ensemble = generateEnsemble(
      molecule,
      100
    );

    BOOST_CHECK(
      TemplateMagic::all_of(
        TemplateMagic::map(
          ensemble,
          [&](const auto& positions) -> bool {
            return molecule.stereocenters
            == molecule.getAdjacencyList().inferStereocentersFromPositions(
              positions
            );
          }
        )
      )
    );
  }
}
