#define BOOST_TEST_MODULE DGDistanceBoundsMatrixTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "WriteMatrix.h"
#include "BoundsFromSymmetry.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/generateConformation.h"

BOOST_AUTO_TEST_CASE( DistanceBoundsTests ) {
  using namespace MoleculeManip;
  using namespace MoleculeManip::DistanceGeometry;

  unsigned N = 4;


  Eigen::MatrixXd testMatrix(N, N);
  testMatrix <<   0,   1, 100,   1,
                  1,   0,   1, 100,
                0.5,   1,   0,   1,
                  1, 0.5,   1,   0;

  DistanceBoundsMatrix testBounds(testMatrix);

  testBounds.smooth();

  BOOST_CHECK(testBounds.upperBound(0, 2) == 2 && testBounds.upperBound(1, 3) == 2);

  auto distancesMatrix = testBounds.generateDistanceMatrix(
    MetrizationOption::off
  );

  for(unsigned i = 0; i < N; i++) {
    for(unsigned j = i + 1; j < N; j++) {
      BOOST_CHECK(
        distancesMatrix(i, j) <= testBounds.upperBound(i, j)
        && distancesMatrix(i, j) >= testBounds.lowerBound(i, j)
      );
    }
  }
}

Eigen::MatrixXd boundsToSlack(const Eigen::MatrixXd& bounds) {
  const unsigned N = bounds.cols();
  Eigen::MatrixXd slack(N, N);
  slack.setZero();

  for(unsigned i = 0; i < N; i++) {
    for(unsigned j = i + 1; j < N; j++) {
      slack(i, j) = bounds(i, j) - bounds(j, i);
    }
  }

  return slack;
}

BOOST_AUTO_TEST_CASE( boundsFromSymmetryTests ) {
  for(const auto& symmetryName : Symmetry::allNames) {
    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    auto molecule = DGDBM::symmetricMolecule(symmetryName);
    auto boundsMatrix = MoleculeManip::DistanceGeometry::gatherDGInformation(molecule).distanceBounds;

    writeMatrix(
      "pre-"s + spaceFreeName,
      boundsToSlack(
        boundsMatrix.access()
      )
    );

    boundsMatrix.smooth();

    writeMatrix(
      "post-"s + spaceFreeName,
      boundsToSlack(
        boundsMatrix.access()
      )
    );
  }
}
