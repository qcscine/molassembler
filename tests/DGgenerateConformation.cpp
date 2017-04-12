#define BOOST_TEST_MODULE DGRefinementProblemTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BoundsFromSymmetry.h"
#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/MetricMatrix.h"
#include "Log.h"

BOOST_AUTO_TEST_CASE( propagator ) {
  Log::level = Log::Level::None;
  Log::particulars = {Log::Particulars::PrototypePropagatorDebugInfo};

  using namespace MoleculeManip::DistanceGeometry;
  using namespace MoleculeManip::Stereocenters;

  auto molecule = DGDBM::symmetricMolecule(Symmetry::Name::SquareAntiPrismatic);
  auto DGInfo = gatherDGInformation(molecule);

  auto distances = DGInfo.distanceBounds.generateDistanceMatrix();

  detail::PrototypePropagator propagator {distances};

  for(const auto& prototype : DGInfo.chiralityConstraintPrototypes) {
    // May not fail through assert failure
    propagator(prototype);
  }
}
