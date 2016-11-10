#ifndef INCLUDE_GRAPH_DISTANCE_MATRIX_H
#define INCLUDE_GRAPH_DISTANCE_MATRIX_H

#include "AdjacencyMatrix.h"
#include <Eigen/Core>

namespace MoleculeManip {

class GraphDistanceMatrix {
private:
  Eigen::MatrixXd _matrix;

  void _transformToDistances();
  void _copyInRow(
    const AtomIndexType& sourceRow,
    const AtomIndexType& targetRow
  );

public:
  const unsigned N;

  GraphDistanceMatrix() = delete;
  GraphDistanceMatrix(const AdjacencyMatrix& adjacencyMatrix);
  GraphDistanceMatrix(AdjacencyMatrix&& adjacencyMatrix);

  std::vector<
    std::vector<AtomIndexType>
  > extractChains(
    const AtomIndexType& i,
    const AtomIndexType& j
  );

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
