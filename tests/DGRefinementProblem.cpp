#define BOOST_TEST_MODULE DGRefinementProblemTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "cppoptlib/meta.h"

#include "BoundsFromSymmetry.h"
#include "IO.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

Eigen::MatrixXd squareMatrixFromRowWiseVector(
  std::vector<double>&& vectorTemporary
) {
  unsigned N = sqrt(vectorTemporary.size());
  if(N*N != vectorTemporary.size()) {
    throw std::logic_error(
      "Asked to create square Matrix from non-square length vector!"
    );
  }

  Eigen::MatrixXd matrix(N, N);

  for(unsigned i = 0; i < N; i++) {
    for(unsigned j = 0; j < N; j++) {
      matrix(i, j) = vectorTemporary[i * N + j];
    }
  }
  
  return matrix;
}

std::vector<Eigen::MatrixXd> testMatrices {
  squareMatrixFromRowWiseVector({
      0,   1,   2,   1,
      1,   0,   1,   2,
    0.5,   1,   0,   1,
      1, 0.5,   1,   0
  })
};

BOOST_AUTO_TEST_CASE( DGRefinementProblemCorrectness ) {

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
