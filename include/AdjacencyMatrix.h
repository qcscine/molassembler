#ifndef INCLUDE_ADJACENCY_MATRIX_H
#define INCLUDE_ADJACENCY_MATRIX_H

#include "AdjacencyList.h"
#include <Eigen/Core>

namespace MoleculeManip {

class AdjacencyMatrix {
private:
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> _matrix;

public:
  const unsigned N;

  AdjacencyMatrix() = delete;
  AdjacencyMatrix(const AdjacencyList& adjacencyList);

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& getMatrixRef() {
    return _matrix;
  }

  bool& operator () (
    const AtomIndexType& i,
    const AtomIndexType& j
  );
  bool operator () (
    const AtomIndexType& i,
    const AtomIndexType& j
  ) const;
};

}

#endif
