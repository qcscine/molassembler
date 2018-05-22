#define BOOST_TEST_MODULE MoleculeSpatialModelTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "BoundsFromSymmetry.h"

BOOST_AUTO_TEST_CASE(dumpDebugInfo) {
  using namespace molassembler::DistanceGeometry;

  for(const auto& symmetryName: Symmetry::allNames) {
    auto molecule = DGDBM::asymmetricMolecule(symmetryName);

    // Default-assign any unassigned stereocenters
    for(const auto& stereocenterPtr : molecule.getStereocenterList()) {
      if(!stereocenterPtr -> assigned()) {
        molecule.assignStereocenter(
          stereocenterPtr->involvedAtoms().front(),
          0
        );
      }
    }

    MoleculeSpatialModel spatialModel {molecule};

    std::cout << Symmetry::name(symmetryName) << std::endl;
    spatialModel.dumpDebugInfo();

    std::cout << "Resulting bounds matrix:" << std::endl;
    DistanceBoundsMatrix bounds {molecule, spatialModel.makeBoundList()};

    bounds.smooth();

    std::cout << bounds.access() << std::endl;
  }
}
