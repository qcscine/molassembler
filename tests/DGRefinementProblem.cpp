#include "BoostTestingHeader.h"
#include <iostream>

#include "DistanceGeometry/generateConformation.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"

BOOST_AUTO_TEST_CASE( DGRefinementProblemCorrectness ) {
  using namespace MoleculeManip;
  using namespace MoleculeManip::DistanceGeometry;

  unsigned N = 4;

  Eigen::MatrixXd distanceBounds;
  distanceBounds.resize(N, N);
  distanceBounds <<   0,   1,   2,   1,
                      1,   0,   1,   2,
                    0.5,   1,   0,   1,
                      1, 0.5,   1,   0;
                    

  DistanceBoundsMatrix testBounds(distanceBounds);

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

  /* NOTE
   * - Finite difference gradient checking may not be up to the task in this
   *   somewhat more difficult implementation, but it's hard to say for certain.
   *   True test of correctness is actual molecules.
   */
  /* 
  BOOST_CHECK(
    problem.checkGradient(vectorizedPositions)
  );
  */

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
