#ifndef INCLUDE_TRIANGULAR_MATRIX_H
#define INCLUDE_TRIANGULAR_MATRIX_H

#include <Eigen/Core>

namespace MoleculeManip {

class TriangularMatrix {
private:
  Eigen::MatrixXd _matrix;

public:
  const unsigned N;

  TriangularMatrix(Eigen::MatrixXd matrix);

  decltype(_matrix(1, 2))& upper(
    const unsigned& i,
    const unsigned& j
  );
  decltype(_matrix(1, 2))& lower(
    const unsigned& i,
    const unsigned& j
  );
  double upper(
    const unsigned& i,
    const unsigned& j
  ) const;
  double lower(
    const unsigned& i,
    const unsigned& j
  ) const;

  Eigen::MatrixXd getMatrix();
};

} // eo namespace MoleculeManip

#endif
