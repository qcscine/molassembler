#include "AdjacencyMatrix.h"

namespace MoleculeManip {

AdjacencyMatrix::AdjacencyMatrix(const AdjacencyList& adjacencyList) 
: N(adjacencyList.size()) {
  _matrix.resize(N, N);
  _matrix.triangularView<Eigen::StrictlyUpper>().setZero();

  for(AtomIndexType i = 0; i < N; i++) {
    for(const auto& adjacentIndex: adjacencyList[i]) {
      if(i < adjacentIndex) {
        this->operator()(i, adjacentIndex) = 1;
      }
    }
  }
}

double& AdjacencyMatrix::operator () (
  const AtomIndexType& i,
  const AtomIndexType& j
) {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}

double AdjacencyMatrix::operator () (
  const AtomIndexType& i,
  const AtomIndexType& j
) const {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}


} // eo namespace
