#define BOOST_TEST_MODULE DGDistanceBoundsMatrixTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "BoundsFromSymmetry.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/generateConformation.h"
#include "temple/Enumerate.h"

BOOST_AUTO_TEST_CASE( DistanceBoundsTests ) {
  using namespace molassembler;
  using namespace molassembler::DistanceGeometry;

  unsigned N = 4;


  Eigen::MatrixXd testMatrix(N, N);
  testMatrix <<   0,   1, 100,   1,
                  1,   0,   1, 100,
                0.5,   1,   0,   1,
                  1, 0.5,   1,   0;

  DistanceBoundsMatrix testBounds(testMatrix);

  testBounds.smooth();

  BOOST_CHECK(testBounds.upperBound(0, 2) == 2 && testBounds.upperBound(1, 3) == 2);

  auto distancesMatrixResult = testBounds.makeDistanceMatrix();
  BOOST_REQUIRE_MESSAGE(distancesMatrixResult, distancesMatrixResult.error().message());
  auto distancesMatrix = distancesMatrixResult.value();

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

std::string flattenMatrix(const Eigen::MatrixXd& matrix) {
  std::string retString;
  
  for(unsigned i = 0; i < matrix.rows(); i++) {
    for(unsigned j = 0; j < matrix.cols(); j++) {
      retString += std::to_string(matrix(i, j));
      if(j != matrix.cols() - 1) retString += ","s;
    }
  }

  return retString;
}

BOOST_AUTO_TEST_CASE( boundsFromSymmetryTests ) {
  std::ofstream outStream("DGDistanceBoundsMatrix-matrix-smoothing.csv");

  for(const auto& enumPair : enumerate(Symmetry::allNames)) {
    const auto& symmetryName = enumPair.value;

    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    auto molecule = DGDBM::symmetricMolecule(symmetryName);
    auto info = molassembler::DistanceGeometry::gatherDGInformation(molecule);

    molassembler::DistanceGeometry::DistanceBoundsMatrix boundsMatrix {
      molecule,
      info.boundList
    };

    outStream << enumPair.index << "," << 0 << ","
      << flattenMatrix(
        boundsToSlack(
          boundsMatrix.access()
        )
      ) << "\n";

    boundsMatrix.smooth();

    outStream << enumPair.index << "," << 1 << ","
      << flattenMatrix(
        boundsToSlack(
          boundsMatrix.access()
        )
      ) << "\n";
  }

  outStream.close();
}
