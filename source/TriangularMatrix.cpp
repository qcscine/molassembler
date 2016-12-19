#include "TriangularMatrix.h"

#include <cassert>

namespace MoleculeManip {

TriangularMatrix::TriangularMatrix(Eigen::MatrixXd matrix) : 
  _matrix(matrix),
  N(matrix.rows()) {
  assert(matrix.rows() == matrix.cols());
}

decltype(TriangularMatrix::_matrix(1, 2))& TriangularMatrix::upper(
  const unsigned& i,
  const unsigned& j
) {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}

decltype(TriangularMatrix::_matrix(1, 2))& TriangularMatrix::lower(
  const unsigned& i,
  const unsigned& j
) {
  return _matrix(
    std::max(i, j),
    std::min(i, j)
  );
}

double TriangularMatrix::upper(
  const unsigned& i,
  const unsigned& j
) const {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}

double TriangularMatrix::lower(
  const unsigned& i,
  const unsigned& j
) const {
  return _matrix(
    std::max(i, j),
    std::min(i, j)
  );
}

Eigen::MatrixXd TriangularMatrix::getMatrix() {
  return _matrix;
}

} // eo namespace
