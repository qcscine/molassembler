#include "BoostTestingHeader.h"
#include <iostream>

#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/MetricMatrix.h"

BOOST_AUTO_TEST_CASE( MetricMatrixTests ) {
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
    
  auto distanceMatrix = testBounds.generateDistanceMatrix(
    MetrizationOption::off
  );

  std::cout << "Distance Matrix: " << std::endl;
  std::cout << distanceMatrix << std::endl << std::endl;

  MetricMatrix metric(
    std::move(distanceMatrix)
  );

  std::cout << "Metric Matrix: " << std::endl;
  std::cout << metric << std::endl << std::endl;

  auto embedded = metric.embed(EmbeddingOption::threeDimensional);

  std::cout << "Embedded positions: " << std::endl;
  std::cout << embedded << std::endl << std::endl;
}
