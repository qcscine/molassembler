#define BOOST_TEST_MODULE DGRefinementProblemTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/MetricMatrix.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "cppoptlib/meta.h"
#include "DistanceGeometry/DGRefinementProblem.h"

#include "BoundsFromSymmetry.h"
#include "IO.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

BOOST_AUTO_TEST_CASE( cppoptlibGradientCorrectness ) {
  // Generate a wide array of DGRefinementProblems and check their gradients

  for(const auto& symmetryName: Symmetry::allNames) {
    auto molecule = DGDBM::symmetricMolecule(symmetryName);
    auto DGInfo = gatherDGInformation(molecule);

    auto distances = DGInfo.distanceBounds.generateDistanceMatrix(
      MetrizationOption::off
    );

    MetricMatrix metric(distances);

    auto embedded = metric.embed(EmbeddingOption::threeDimensional);

    Eigen::VectorXd vectorizedPositions(
      Eigen::Map<Eigen::VectorXd>(
        embedded.data(),
        embedded.cols() * embedded.rows()
      )
    );

    auto chiralityConstraints = TemplateMagic::map(
      DGInfo.chiralityConstraintPrototypes,
      detail::PrototypePropagator {distances}
    );

    DGRefinementProblem<double> problem(
      chiralityConstraints,
      DGInfo.distanceBounds
    );

    BOOST_CHECK(
      problem.checkGradient(vectorizedPositions)
    );
  }
}
