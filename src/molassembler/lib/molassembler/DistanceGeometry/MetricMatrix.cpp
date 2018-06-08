#include "DistanceGeometry/MetricMatrix.h"
#include "detail/SharedTypes.h"

#include <Eigen/Eigenvalues>

/* TODO
 * - Is numerical accuracy maybe an issue somewhere here?
 */

// C++17 alter to namespace molassembler::DistanceGeometry {
namespace molassembler {

namespace DistanceGeometry {

void MetricMatrix::_constructFromTemporary(Eigen::MatrixXd&& distances) {
  // Beware, only strict upper triangle of distances contains anything

  const AtomIndexType N = distances.rows();

  // resize and null-initialize matrix
  _matrix.resize(N, N);
  _matrix.setZero();

  // Declare D0 vector
  Eigen::VectorXd D0(N);
  /* D_{i}² =   (1/N) * sum_{j}(distances[i, j]²)
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
    D0[i] = firstTerm - doubleSumTerm;
  }

  /* The diagonal of G is just D0!
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
   * to embed the metric matrix
   */
  for(AtomIndexType i = 0; i < N; i++) {
    for(AtomIndexType j = i + 1; j < N; j++) {
      _matrix(j, i) = (
        D0[i] + D0[j] - distances(i, j)
      ) / 2.0;
    }
  }
}

MetricMatrix::MetricMatrix(Eigen::MatrixXd&& matrix) {
  _constructFromTemporary(std::forward<Eigen::MatrixXd>(matrix));
}

MetricMatrix::MetricMatrix(const Eigen::MatrixXd& matrix) {
  auto copy = matrix;
  _constructFromTemporary(std::move(copy));
}

const Eigen::MatrixXd& MetricMatrix::access() const {
  return _matrix;
}

Eigen::MatrixXd MetricMatrix::embed() const {
  const unsigned dimensionality = 4;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(_matrix);

  // reverse because smallest are listed first by Eigen
  Eigen::VectorXd eigenValues = eigenSolver.eigenvalues().reverse();
  auto numEigenValues = eigenValues.size();

  eigenValues.conservativeResize(dimensionality); // dimensionality x 1

  // If we have upscaled, the new elements must be zero!
  if(numEigenValues < dimensionality) {
    for(unsigned i = numEigenValues; i < dimensionality; i++) {
      eigenValues(i) = 0;
    }
  }

  // If any eigenvalues in the vector are negative, set them to 0
  for(unsigned i = 0; i < dimensionality; i++) {
    if(eigenValues(i) < 0) {
      eigenValues(i) = 0;
    }
  }

  // take square root of eigenvalues
  eigenValues = eigenValues.cwiseSqrt();

  Eigen::MatrixXd L;
  L.resize(dimensionality, dimensionality);
  L.setZero();
  L.diagonal() = eigenValues;

  Eigen::MatrixXd V = eigenSolver.eigenvectors();
  // Eigen has its own concept of rows and columns, I would have thought it's
  // columns. But tests have shown it has to be row-wise.
  V.rowwise().reverseInPlace();
  V.conservativeResize(
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

bool MetricMatrix::operator == (const MetricMatrix& other) const {
  return _matrix == other._matrix;
}

std::ostream& operator << (
  std::ostream& os,
  const MetricMatrix& metricMatrix
) {
  os << metricMatrix._matrix;
  return os;
}

} // namespace DistanceGeometry

} // namespace molassembler
