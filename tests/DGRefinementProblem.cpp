#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DGRefinementProblemTests
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "cppoptlib/meta.h"

#include "BoundsFromSymmetry.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

void writeFile(
  const bool& optimized,
  const std::string symmetryString,
  const unsigned& structNum,
  const Eigen::VectorXd& vectorizedPositions
) {
  unsigned N = vectorizedPositions.size() / 3;
  assert(3 * N == vectorizedPositions.size());

  std::vector<std::string> elementNames = {
    "Ru",
    "H",
    "F",
    "Cl",
    "Br",
    "I",
    "N",
    "C",
    "O",
    "S",
    "P"
  };

  std::string filename;

  if(optimized) filename = "opt-"s;
  else filename = "gen-"s;

  filename += symmetryString + std::to_string(structNum) + ".xyz"s;

  std::ofstream outStream(filename.c_str());

  outStream << std::setprecision(7) << std::fixed;

  outStream << N << std::endl
    << "Energy = " << std::endl;

  for(unsigned i = 0; i < N; i++) {
    if(elementNames[i].size() == 1) outStream << elementNames[i] << " ";
    else outStream << elementNames[i];

    outStream << std::setw(13) << vectorizedPositions[3 * i + 0];
    outStream << std::setw(13) << vectorizedPositions[3 * i + 1];
    outStream << std::setw(13) << vectorizedPositions[3 * i + 2];

    if(i != N - 1) outStream << std::endl;
  }

  outStream.close();
}

void writeErrorValues(
  const std::string& symmetryString,
  const std::vector<double> errorValues
) {
  std::ofstream outStream(symmetryString+".csv"s);

  for(const auto& value: errorValues) {
    outStream << value << std::endl;
  }

  outStream.close();
}


/* NOTES to current state of optimized output structures
 *
 * - Initial problem is gone, was an issue with the generation of the test
 *   matrices
 * - Found one bug in the refinement problem (a missing square), heavily
 *   decreases optimization time
 * - Found another bug that caused generated structures to be considerably
 *   expanded, although well reflective of the overall symmetries (a missing
 *   sqrt in the embedding procedure).
 * - Now high-symmetry geometries are maybe somewhat over-constrained, leading 
 *   to messy structures. Avenues to try:
 *
 *   - Add more variance to 1-3 distance bounds
 *     -> did nothing
 *
 *   - Add more sanity tests to the various symmetries' angle functions
 *     -> Found mistakes in the angle functions of pent-bipy and sq-antiprism, 
 *        fixing the issue
 *        
 * - Refinement cannot reach a global minimum in almost all cases. Probably due
 *   to triangle inequality bound violations. Fix by introducing metrization.
 */

BOOST_AUTO_TEST_CASE(distanceBoundsGeneratedArePlausible) {
  // Make a problem solver
  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGSolver;

  // Set stop criteria
  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  //stopCriteria.iterations = 1000;
  stopCriteria.fDelta = 1e-5;
  DGSolver.setStopCriteria(stopCriteria);

  const unsigned nStructures = 100;

  for(const auto& symmetryName : Symmetry::allNames) {
    // Make a space-free string from the name
    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    // Generate distance bounds
    auto distanceBoundsMatrix = DGDBM::distanceBoundsFromSymmetry(
      symmetryName,
      DGDBM::DistancesOption::Uniform
    );

    /*std::cout << "Sample distances matrix for symmetry '" 
      << Symmetry::name(symmetryName) << std::endl
      << distanceBoundsMatrix.generateDistanceMatrix(
        MetrizationOption::off
      ) << std::endl;*/

    std::vector<double> refinedErrorValues;

    for(unsigned structNum = 0; structNum < nStructures; structNum++) {

      // Calculate metric matrix from selected distances
      MetricMatrix metricMatrix(
        distanceBoundsMatrix.generateDistanceMatrix(
          MetrizationOption::off
        )
      );

      // Embed
      auto embeddedPositions = metricMatrix.embed(EmbeddingOption::threeDimensional);

      // Vectorize
      Eigen::VectorXd vectorizedPositions(
        Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      );

      // Write the unoptimized result to a file
      writeFile(
        false,
        spaceFreeName,
        structNum,
        vectorizedPositions
      );

      // Create the RefinementProblem
      DGRefinementProblem<double> problem(
        std::vector<ChiralityConstraint>({}), // no chirality constraints
        distanceBoundsMatrix
      );

      // Run the minimization
      DGSolver.minimize(problem, vectorizedPositions);

      // Save end value in Problem
      refinedErrorValues.emplace_back(
        problem.value(
          vectorizedPositions
        )
      );

      // Write the result to a file
      writeFile(
        true,
        spaceFreeName,
        structNum,
        vectorizedPositions
      );
    }

    // Write refined error values to file
    writeErrorValues(
      spaceFreeName,
      refinedErrorValues
    );
  }
}

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
