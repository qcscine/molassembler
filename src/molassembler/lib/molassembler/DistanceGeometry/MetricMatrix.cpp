/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/DistanceGeometry/MetricMatrix.h"
#include "molassembler/Types.h"

#include <Eigen/Eigenvalues>

namespace Scine {

namespace molassembler {

namespace distance_geometry {

void MetricMatrix::_constructFromTemporary(Eigen::MatrixXd&& distances) {
  /* We have to be a little careful since only strict upper triangle of
   * distances contains anything of use to us.
   *
   * So we have to make sure the first index is always smaller than the second
   * to reference the correct entry.
   */

  const AtomIndex N = distances.rows();

  /* Resize our underlying matrix. There is no need to zero-initialize since all
   * parts of the lower triangle (including the diagonal) are overwritten.
   */
  _matrix.resize(N, N);

  /* We need to accomplish the following:
   *
   * Every entry in the lower triangle of _matrix (including the diagonal)
   * needs to be set according to the result of the following equations:
   *
   * D0[i]² =   (1/N) * sum_{j}(distances[i, j]²)
   *          - (1/(N²)) * sum_{j < k}(distances[j, k]²)
   *
   *  (The second term is independent of i and can be precalculated!)
   *
   * _matrix[i, j] = ( D0[i]² + D0[j]² - distances[i, j]² ) / 2
   *
   * On the diagonal, where i == j:
   * _matrix[i, i] = ( D0[i]² + D0[i]² - distances[i, i]² ) / 2
   *                     ^--------^      ^-------------^
   *                       equal               =0
   *
   * -> _matrix[i, i] = D0[i]²
   *
   * So, we can store all of D0 immediately on _matrix's diagonal and perform
   * the remaining transformation afterwards.
   */

  // Since we need squares EVERYWHERE, just square the whole distances matrix
  distances = distances.cwiseProduct(distances);

  double doubleSumTerm = 0;
  for(AtomIndex j = 0; j < N; ++j) {
    for(AtomIndex k = j + 1; k < N; ++k) {
      doubleSumTerm += distances(j, k);
    }
  }
  doubleSumTerm /= N * N;

  for(AtomIndex i = 0; i < N; ++i) {
    // compute first term
    double firstTerm = 0;
    for(AtomIndex j = 0; j < N; ++j) {
      if(i == j) {
        continue;
      }

      firstTerm += distances(
        std::min(i, j),
        std::max(i, j)
      );
    }
    firstTerm /= N;

    // assign as difference, no need to sqrt, we need the squares in a moment
    _matrix.diagonal()(i) = firstTerm - doubleSumTerm;
  }

  /* Write off-diagonal elements into lower triangle
   * Why the lower triangle? Because that is the only part of the matrix
   * referenced by Eigen's SelfAdjointEigenSolver, which we will use in a bit
   * to embed the metric matrix
   */
  for(AtomIndex i = 0; i < N; i++) {
    for(AtomIndex j = i + 1; j < N; j++) {
      _matrix(j, i) = (
        // D0[i]²             + D0[j]²                - d(i, j)²
        _matrix.diagonal()(i) + _matrix.diagonal()(j) - distances(i, j)
      ) / 2.0;
    }
  }
}

MetricMatrix::MetricMatrix(Eigen::MatrixXd distanceMatrix) {
  _constructFromTemporary(std::move(distanceMatrix));
}

const Eigen::MatrixXd& MetricMatrix::access() const {
  return _matrix;
}

Eigen::MatrixXd MetricMatrix::embed() const {
  return embedWithFullDiagonalization();
}

Eigen::MatrixXd MetricMatrix::embedWithFullDiagonalization() const {
  constexpr unsigned dimensionality = 4;

  // SelfAdjointEigenSolver only references the lower triangle
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(_matrix);

  Eigen::VectorXd eigenvalues = eigenSolver.eigenvalues();

  // Construct L
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(dimensionality, dimensionality);
  // We want the algebraically largest eigenvalues (up to four, if present)
  unsigned numEigenvalues = std::min(
    static_cast<unsigned>(eigenvalues.size()),
    dimensionality
  );
  for(unsigned i = 0; i < numEigenvalues; ++i) {
    /* Since Eigen stores them in increasing order, we have to fetch the
     * algebraically largest from the back. We only want to use the eigenpair
     * if the eigenvalue is greater than zero.
     */
    if(eigenvalues(eigenvalues.size() - i - 1) > 0) {
      L.diagonal()(i) = std::sqrt(
        eigenvalues(eigenvalues.size() - i - 1)
      );
    }
  }

  // V is initially N x N
  Eigen::MatrixXd V = eigenSolver.eigenvectors();
  /* Again, eigenvectors are sorted in increasing corresponding eigenvalues
   * algebraic value. We have to reverse them to match the ordering in L.
   */
  V.rowwise().reverseInPlace();
  V.conservativeResize(V.rows(), dimensionality);
  // now N x dimensionality

  /* Calculate X = VL
   * (N x 4) · (4 x 4) -> (N x 4), but we want (4 x N), so we transpose
   */
  return (V * L).transpose();
}

bool MetricMatrix::operator == (const MetricMatrix& other) const {
  return _matrix == other._matrix;
}

} // namespace distance_geometry

} // namespace molassembler

} // namespace Scine
