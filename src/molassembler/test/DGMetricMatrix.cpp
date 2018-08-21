#define BOOST_TEST_MODULE DGMetricMatrixTests
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"
#include "chemical_symmetries/Symmetries.h"
#include "Eigen/Eigenvalues"
#include "temple/Stringify.h"

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/MetricMatrix.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"
#include "molassembler/detail/StdlibTypeAlgorithms.h"
#include "molassembler/IO.h"
#include "molassembler/Options.h"

#include <iostream>
#include <random>

using namespace molassembler;
using namespace molassembler::DistanceGeometry;

Eigen::MatrixXd reorder(
  const Eigen::MatrixXd& sourceMatrix,
  const std::vector<unsigned>& reorderSequence
) {
  Eigen::MatrixXd retMatrix(reorderSequence.size(), reorderSequence.size());

  /* e.g. reorderSequence is {4, 2, 1, 3}
   *  and input is
   *
   *     1  2  3  4          4  2  1  3
   *
   * 1   -  1  2  3      4   -  5  3  6
   * 2   -  -  4  5  ->  2   -  -  1  4
   * 3   -  -  -  6  ->  1   -  -  -  2
   * 4   -  -  -  -      3   -  -  -  -
   *
   * reverse:
   *
   *     1  2  3  4          3  2  4  1
   *
   * 1   -  5  3  6      3   -  1  2  3
   * 2   -  -  1  4  ->  2   -  -  4  5
   * 3   -  -  -  2  ->  4   -  -  -  6
   * 4   -  -  -  -      1   -  -  -  -
   */


  for(unsigned i = 0; i < reorderSequence.size(); i++) {
    // Diagonal
    retMatrix(i, i) = sourceMatrix(
      reorderSequence[i],
      reorderSequence[i]
    );

    // Triangles
    for(unsigned j = i + 1; j < reorderSequence.size(); j++) {
      // Upper
      retMatrix(i, j) = sourceMatrix(
        std::min(
          reorderSequence[i],
          reorderSequence[j]
        ),
        std::max(
          reorderSequence[i],
          reorderSequence[j]
        )
      );

      // Lower
      retMatrix(j, i) = sourceMatrix(
        std::max(
          reorderSequence[i],
          reorderSequence[j]
        ),
        std::min(
          reorderSequence[i],
          reorderSequence[j]
        )
      );
    }
  }

  return retMatrix;
}

std::vector<unsigned> inverseReorderSequence(
  const std::vector<unsigned>& reorderSequence
) {
  /* reorder:
   *
   * 1 -> 4
   * 2 -> 2
   * 3 -> 1
   * 4 -> 3
   *
   * inverse:
   *
   * 1 -> 3 (what position 1 was at in the reorderSequence)
   * 2 -> 2 (...)
   * 3 -> 4
   * 4 -> 1
   */
  std::vector<unsigned> returnSequence (reorderSequence.size());

  for(unsigned i = 0; i < reorderSequence.size(); i++) {
    auto findIter = std::find(
      reorderSequence.begin(),
      reorderSequence.end(),
      i
    );
    if(findIter == reorderSequence.end()) {
      throw std::logic_error("Could not find i in reorderSequence!");
    }

    returnSequence[i] = findIter - reorderSequence.begin();
  }

  return returnSequence;
}

std::vector<unsigned> randomReorderingSequence(const unsigned length) {
  std::vector<unsigned> reorderSequence (length);

  std::iota(
    reorderSequence.begin(),
    reorderSequence.end(),
    0
  );

  prng.shuffle(reorderSequence);

  return reorderSequence;
}

BOOST_AUTO_TEST_CASE( reorderingWorks ) {
  bool allPassed = true;
  const unsigned N = 10;

  for(unsigned i = 0; i < 100; i++) {
    const Eigen::MatrixXd testMatrix = Eigen::MatrixXd::Random(N, N);

    auto reorderingOrder = randomReorderingSequence(N);

    auto result = reorder(
      reorder(
        testMatrix,
        reorderingOrder
      ),
      inverseReorderSequence(
        reorderingOrder
      )
    );

    if(testMatrix != result) {
      allPassed=false;

      std::cout << "Reordering reversibility failed! Original:\n"
        << testMatrix << std::endl
        << "Reordering sequence: " << temple::stringify(reorderingOrder) << std::endl
        << "Unreordering sequence: " << temple::stringify(inverseReorderSequence(reorderingOrder))
        << "Computed result: " << result << std::endl;

      break;
    }
  }

  BOOST_CHECK(allPassed);
}

void showEmbedding(const MetricMatrix& metricMatrix) {
  auto dimensionality = 4u;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(metricMatrix.access());

  // reverse because smallest are listed first by Eigen
  Eigen::VectorXd eigenValues = eigenSolver.eigenvalues().reverse();

  auto numEigenValues = eigenValues.size();

  std::cout << "Reversed eigenvalues of metric matrix:\n" << eigenValues << std::endl;

  eigenValues.conservativeResize(dimensionality); // dimensionality x 1

  std::cout << "reversed and resized eigenvalues of metric matrix:\n" << eigenValues << std::endl;

  // If we have upscaled the eigenvalues, any new elements must be set to zero
  if(numEigenValues < dimensionality) {
    for(unsigned i = numEigenValues; i < dimensionality; i++) {
      eigenValues(i) = 0;
    }

    std::cout << "We had fewer eigenvalues than dimensionality requires, so the new elements must be set zero:\n"
      << eigenValues << std::endl;
  }

  // If any eigenvalues in the vector are negative, set them to 0
  for(unsigned i = 0; i < dimensionality; i++) {
    if(eigenValues(i) < 0) {
      eigenValues(i) = 0;
    }
  }

  std::cout << "After setting any negative eigenvalues to zero:\n" << eigenValues << std::endl;

  // take square root of eigenvalues
  eigenValues = eigenValues.cwiseSqrt();

  std::cout << "sqrt-ed eigenvales of metric matrix: \n" << eigenValues << std::endl;

  Eigen::MatrixXd L;
  L.resize(dimensionality, dimensionality);
  L.setZero();
  L.diagonal() = eigenValues;

  std::cout << "Resulting L matrix:\n" << L << std::endl;

  Eigen::MatrixXd V = eigenSolver.eigenvectors();

  std::cout << "Eigenvectors:\n" << V << std::endl;
  // Eigen has its own concept of rows and columns, I would have thought it's
  // columns. But tests have shown it has to be row-wise.
  V.rowwise().reverseInPlace();

  std::cout << "Post row-wise reverse:\n" << V << std::endl;

  V.conservativeResize(
    V.rows(),
    dimensionality
  ); // now Natoms x dimensionality

  std::cout << "Post reduce to dimensionality:\n" << V << std::endl;

  /* V * L
   * (Natoms x dimensionality) · (dimensionality x dimensionality )
   * -> (Natoms x dimensionality)
   * transpose (V * L)
   * -> dimensionality x Natoms
   */
  std::cout << "Resulting positions matrix:\n" << (V * L).transpose() << std::endl;
}

BOOST_AUTO_TEST_CASE( constructionIsInvariantUnderOrderingSwap ) {
  boost::filesystem::path filesPath("test_files/ez_stereocenters");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;

    auto molecule = IO::read(currentFilePath.string());

    auto DGData = DistanceGeometry::gatherDGInformation(molecule);

    DistanceBoundsMatrix distanceBounds {
      molecule,
      DGData.bounds
    };

    // choose a random reordering
    auto reorderSequence = randomReorderingSequence(
      molecule.numAtoms() + 1
    );

    // get a distances matrix from the bounds
    auto distancesMatrixResult = distanceBounds.makeDistanceMatrix();
    if(!distancesMatrixResult) {
      BOOST_FAIL(distancesMatrixResult.error().message());
    }
    auto distancesMatrix = distancesMatrixResult.value();

    // reorder the distance Bounds matrix
    auto reorderedDM = reorder(distancesMatrix, reorderSequence);

    // generate a metric matrix for both
    auto originalMetric = MetricMatrix(std::move(distancesMatrix));
    auto reorderedMetric = MetricMatrix(std::move(reorderedDM));

    /* Metric Matrix of original distances has to be identical to
     * un-reordered metric matrix of reordered distances
     */
    auto revert = reorder(
      reorderedMetric.access(),
      inverseReorderSequence(reorderSequence)
    );

    if(
      !originalMetric.access().isApprox(
        revert,
        1e-7
      )
    ) {
      std::cout << "Failed reordering test for "
        << currentFilePath.string() << ":" << std::endl
        << "Metric matrix from original distances matrix: " << std::endl
        << originalMetric << std::endl
        << "un-reordered Metric matrix from reordered:"
        << std::endl << revert << std::endl;
      BOOST_FAIL("Reordering test fails!");
      break;
    }

    /*if(symmetryName == Symmetry::Name::SquareAntiPrismatic) {
      showEmbedding(originalMetric);
    }*/

    /*if(nTests == 0) { // once per symmetry
      showEmbedding(originalMetric);
    }*/
  }
}

BOOST_AUTO_TEST_CASE( explicitFromLecture ) {
  // From Algorithms & Programming in C++ lecture
  Eigen::MatrixXd exactDistanceMatrix (4, 4);
  // No need to enter lower triangle
  exactDistanceMatrix << 0, 1, sqrt(2),       1,
                         0, 0,       1, sqrt(2),
                         0, 0,       0,       1,
                         0, 0,       0,       0;

  auto metric = MetricMatrix(exactDistanceMatrix);

  Eigen::MatrixXd expectedMetricMatrix (4, 4);
  expectedMetricMatrix <<  0.5,    0,    0,    0,
                             0,  0.5,    0,    0,
                          -0.5,    0,  0.5,    0,
                             0, -0.5,    0,  0.5;

  BOOST_CHECK(
    metric.access().isApprox(
      expectedMetricMatrix,
      1e-7
    )
  );
}

BOOST_AUTO_TEST_CASE( highSymmetryFailures ) {
  Eigen::MatrixXd sampleLinearDistances(3, 3);
  sampleLinearDistances <<     0, 1.882, 2.196,
                           1.782,     0, 4.078,
                           2.096, 3.878,     0;

  auto metric = MetricMatrix(sampleLinearDistances);

  showEmbedding(metric);
}
