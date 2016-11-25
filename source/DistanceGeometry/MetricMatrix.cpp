#include "DistanceGeometry/MetricMatrix.h"
#include "common_typedefs.h"

#include <Eigen/Eigenvalues>

// C++17 alter to namespace MoleculeManip::DistanceGeometry {
namespace MoleculeManip {

namespace DistanceGeometry {

MetricMatrix::MetricMatrix(Eigen::MatrixXd&& distances) {
  const AtomIndexType N = distances.rows();
  _matrix.resize(N, N);

  Eigen::VectorXd D0(N);
  /* D_{i}² =   (1/N) * sum_{j}(distances[i, j]) 
   *          - (1/(N²)) * sum_{j < k}(distances[j, k]²)
   *
   * The second term is independent of i, precalculate!
   */

  // Since we need squares EVERYWHERE, just square the whole distances matrix
  distances = distances.cwiseProduct(distances);

  double doubleSumTerm = 0;
  for(AtomIndexType j = 0; j < N; j++) {
    for(AtomIndexType k = j + 1; k < N; k++) {
      doubleSumTerm += distances(j, k);
    }
  }
  doubleSumTerm /= N * N;

  for(AtomIndexType i = 0; i < N; i++) {
    // compute first term
    double firstTerm = 0;
    for(AtomIndexType j = 0; j < N; j++) {
      if(i == j) continue;
      firstTerm += distances(
        std::min(i, j),
        std::max(i, j)
      );
    }
    firstTerm /= N;

    // assign as difference, no need to sqrt, we need the squares in a moment
    D0[i] = firstTerm - doubleSumTerm;
  }

  /* Diagonal is just D0! 
   * g[i, j] = ( D0[i]² + D0[j]² - distances[i, j]² ) / 2
   * if i == j:    ^--------^      ^-------------^
   *                 equal               =0
   *
   * -> g[i, i] = D0[i]²
   *
   * since we store the squares, no need to do anything but assign the entire
   * diagonal as the contents of the previous vector. Could also optimize out 
   * the D0 vector, it's actually unnecessary (just store it on the diagonal 
   * immediately).
   */
  _matrix.diagonal() = D0;

  /* Write off-diagonal elements into lower triangle
   * Why the lower triangle? Because that is the only part of the matrix 
   * referenced by Eigen's SelfAdjointEigenSolver, which we will use in a bit
   */
  for(AtomIndexType i = 0; i < N; i++) {
    for(AtomIndexType j = i + 1; j < N; j++) {
      _matrix(j, i) = (
        D0[i] + D0[j] - distances(i, j)
      ) / 2.0;
    }
  }
}

Eigen::MatrixXd MetricMatrix::embed(
  const EmbeddingOption& embedding 
) const {
  unsigned dimensionality = static_cast<unsigned>(embedding);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(_matrix);

  // reverse because smallest are listed first by Eigen
  Eigen::VectorXd eigenValues = eigenSolver.eigenvalues().reverse();
  eigenValues.conservativeResize(dimensionality); // dimensionality x 1

  Eigen::MatrixXd L;
  L.resize(dimensionality, dimensionality);
  L.setZero();
  L.diagonal() = eigenValues;

  Eigen::MatrixXd V = eigenSolver.eigenvectors();
  // Eigen has its own concept of rows and columns, I would have thought it's 
  // columns. But tests have shown it has to be rowwise.
  V = V.rowwise().reverse();
  V.resize(
    V.rows(),
    dimensionality
  ); // now Natoms x dimensionality

  /* V * L
   * (Natoms x dimensionality) · (dimensionality x dimensionality )
   * -> (Natoms x dimensionality)
   * transpose (V * L)
   * -> dimensionality x Natoms
   */
  return (V * L).transpose(); 
}

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip
