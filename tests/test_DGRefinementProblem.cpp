#include <iostream>

#include "DistanceGeometry/generateConformation.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;


BOOST_AUTO_TEST_CASE( DGRefinementProblemCorrectness ) {
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

  MetricMatrix metric(
    testBounds.generateDistanceMatrix(
      MetrizationOption::off
    )
  );

  auto embedded = metric.embed(EmbeddingOption::threeDimensional);

  std::cout << "Embedded positions: " << std::endl;
  std::cout << embedded << std::endl << std::endl;

  Eigen::VectorXd vectorizedPositions(
    Eigen::Map<Eigen::VectorXd>(
      embedded.data(),
      embedded.cols() * embedded.rows()
    )
  );

  std::cout << "Vectorized positions: " << std::endl;
  std::cout << vectorizedPositions << std::endl << std::endl;
  
  std::vector<ChiralityConstraint> constraints;

  DGRefinementProblem<double> problem(
    constraints,
    testBounds
  );

  BOOST_CHECK(
    problem.checkGradient(vectorizedPositions)
  );

  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGConjugatedGradientDescentSolver;

  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  stopCriteria.iterations = 15;
  stopCriteria.fDelta = 1e-5;

  DGConjugatedGradientDescentSolver.setStopCriteria(stopCriteria);
  DGConjugatedGradientDescentSolver.minimize(problem, vectorizedPositions);

  std::cout << "Vectorized positions post minimization: " << std::endl;
  std::cout << vectorizedPositions << std::endl;
  
}
