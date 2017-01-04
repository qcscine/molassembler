#include "BoostTestingHeader.h"
#include <iostream>

#include "DistanceGeometry/DistanceBoundsMatrix.h"

BOOST_AUTO_TEST_CASE( DistanceBoundsTests ) {
  using namespace MoleculeManip;
  using namespace MoleculeManip::DistanceGeometry;

  unsigned N = 4;

  DistanceBoundsMatrix testBounds(N);
  testBounds.upperBound(0, 1) = 1;
  testBounds.upperBound(0, 2) = 2;
  testBounds.upperBound(0, 3) = 1;
  testBounds.upperBound(1, 2) = 1;
  testBounds.upperBound(1, 3) = 2;
  testBounds.upperBound(2, 3) = 1;

  testBounds.lowerBound(0, 1) = 1;
  testBounds.lowerBound(0, 2) = 0.5;
  testBounds.lowerBound(0, 3) = 1;
  testBounds.lowerBound(1, 2) = 1;
  testBounds.lowerBound(1, 3) = 0.5;
  testBounds.lowerBound(2, 3) = 1;

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
