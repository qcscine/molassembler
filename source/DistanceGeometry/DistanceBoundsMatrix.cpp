#include "DistanceGeometry/DistanceBoundsMatrix.h"

#include "StdlibTypeAlgorithms.h"

using namespace MoleculeManip::DistanceGeometry;

void DistanceBoundsMatrix::_setUpperBound(
  const AtomIndexType& i,
  const AtomIndexType& j,
  const double& value
) {
  assert(i != j);
  _matrix(
    std::min(i, j),
    std::max(i, j)
  ) = value;
}

void DistanceBoundsMatrix::_setLowerBound(
  const AtomIndexType& i,
  const AtomIndexType& j,
  const double& value
) {
  assert(i != j);
  _matrix(
    std::max(i, j),
    std::min(i, j)
  ) = value;
}

DistanceBoundsMatrix::DistanceBoundsMatrix(
  const MoleculeManip::Molecule& mol
) {
  // resize the matrix appropriately and set zero
  auto nAtoms = mol.getNumAtoms();
  _matrix.resize(nAtoms, nAtoms);
  _matrix.setZero();

  // populate the matrix with distance constraints
  // initially those from the stereocenters

  // TODO implement constraint getter in Molecule and cache the results!

}
