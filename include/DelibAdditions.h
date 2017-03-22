#ifndef INCLUDE_DELIB_ADDITIONS_H
#define INCLUDE_DELIB_ADDITIONS_H

#include "Types/PositionCollection.h"
#include "common_typedefs.h"

namespace MoleculeManip {

namespace DelibAdditions {

double getAngle(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k
) {
  auto a = positions[i].asEigenVector() - positions[j].asEigenVector(),
       b = positions[k].asEigenVector() - positions[j].asEigenVector();

  return std::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

double getDistance(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j
) {
  return (
    positions[i].asEigenVector()
    - positions[j].asEigenVector()
  ).norm();
}

double toRadians(const double& inDegrees) {
  return M_PI * inDegrees / 180;
}

} // eo namespace DelibAdditions

} // eo namespace MoleculeManip

#endif
