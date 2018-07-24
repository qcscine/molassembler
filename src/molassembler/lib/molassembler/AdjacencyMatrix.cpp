#include "AdjacencyMatrix.h"

#include "Molecule.h"

namespace molassembler {

AdjacencyMatrix::AdjacencyMatrix(const Molecule& molecule)
: N(molecule.numAtoms()) {
  _matrix.resize(N, N);
  _matrix.triangularView<Eigen::StrictlyUpper>().setConstant(false);

  for(AtomIndexType i = 0; i < N; i++) {
    for(const auto& adjacentIndex: molecule[i]) {
      if(i < adjacentIndex) {
        this->operator()(i, adjacentIndex) = true;
      }
    }
  }
}

bool& AdjacencyMatrix::operator () (
  const AtomIndexType i,
  const AtomIndexType j
) {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}

bool AdjacencyMatrix::operator () (
  const AtomIndexType i,
  const AtomIndexType j
) const {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}


} // namespace molassembler
