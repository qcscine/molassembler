#ifndef INCLUDE_ADJACENCY_MATRIX_H
#define INCLUDE_ADJACENCY_MATRIX_H

#include "AdjacencyList.h"
#include <Eigen/Core>

namespace MoleculeManip {

class AdjacencyMatrix {
private:
  Eigen::MatrixXd _matrix;

public:
  const unsigned N;

  AdjacencyMatrix() = delete;
  AdjacencyMatrix(const AdjacencyList& adjacencyList);

  Eigen::MatrixXd& getMatrixRef() {
    return _matrix;
  }

  double& operator () (
    const AtomIndexType& i,
    const AtomIndexType& j
  );
  double operator () (
    const AtomIndexType& i,
    const AtomIndexType& j
  ) const;
};

}

#endif
