#include "DelibHelpers.h"

#include <Eigen/Geometry>

namespace molassembler {

namespace DelibHelpers {

double getDistance(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j
) {
  assert(i < positions.size() && j < positions.size());

  return (
    positions[i].asEigenVector()
    - positions[j].asEigenVector()
  ).norm();
}

double getAngle(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k
) {
  assert(
    i < positions.size()
    && j < positions.size()
    && k < positions.size()
  );

  auto a = positions[i].asEigenVector() - positions[j].asEigenVector(),
       b = positions[k].asEigenVector() - positions[j].asEigenVector();

  return std::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

double getDihedral(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
) {
  assert(
    i < positions.size()
    && j < positions.size()
    && k < positions.size()
    && l < positions.size()
  );

  Eigen::Vector3d a = positions[j].asEigenVector() - positions[i].asEigenVector();
  Eigen::Vector3d b = positions[k].asEigenVector() - positions[j].asEigenVector();
  Eigen::Vector3d c = positions[l].asEigenVector() - positions[k].asEigenVector();

  return std::atan2(
    (
      a.cross(b)
    ).cross(
      b.cross(c)
    ).dot(
      b.normalized()
    ),
    (
      a.cross(b)
    ).dot(
      b.cross(c)
    )
  );
}

double getDihedral(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndexType, 4>& indices
) {
  return getDihedral(
    positions,
    indices[0],
    indices[1],
    indices[2],
    indices[3]
  );
}


double getSignedVolume(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
) {
  assert(
    i < positions.size()
    && j < positions.size()
    && k < positions.size()
    && l < positions.size()
  );

  return (
    positions[i].asEigenVector()
    - positions[l].asEigenVector()
  ).dot(
    (
      positions[j].asEigenVector()
      - positions[l].asEigenVector()
    ).cross(
      positions[k].asEigenVector()
      - positions[l].asEigenVector()
    )
  );
}

double getSignedVolume(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndexType, 4>& indices
) {
  return getSignedVolume(
    positions,
    indices[0],
    indices[1],
    indices[2],
    indices[3]
  );
}

} // namespace DelibHelpers

} // namespace molassembler
