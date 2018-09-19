// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Detail/DelibHelpers.h"

#include <Eigen/Geometry>
#include "temple/Functional.h"

#include <array>

namespace molassembler {

namespace DelibHelpers {

bool validPositionIndices(
  const Delib::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
) {
  /* Delib casts its' underlying vector<Vector3>::size_type to int, so
   * we have to get at the unsigned std::size_t differently. Since vector
   * iterators are RandomAccessIterators, we can get a std::ptrdiff_t in O(1),
   * which is typically bigger than int in most data models.
   */
  auto signedPtrDiffDistance = std::distance(
    std::begin(positions),
    std::end(positions)
  );

  // Casting that to size_t gives us the maximum addressable space back
  std::size_t unsignedPositionCollectionSize = signedPtrDiffDistance;

  return temple::all_of(
    indices,
    [&unsignedPositionCollectionSize](const AtomIndex index) -> bool {
      return index < unsignedPositionCollectionSize;
    }
  );
}

double getDistance(
  const Delib::PositionCollection& positions,
  AtomIndex i,
  AtomIndex j
) {
  assert(
    validPositionIndices(positions, {i, j})
  );

  return (
    positions[i].asEigenVector()
    - positions[j].asEigenVector()
  ).norm();
}

double getAngle(
  const Delib::PositionCollection& positions,
  const AtomIndex i,
  const AtomIndex j,
  const AtomIndex k
) {
  assert(
    validPositionIndices(positions, {i, j, k})
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
  const AtomIndex i,
  const AtomIndex j,
  const AtomIndex k,
  const AtomIndex l
) {
  assert(
    validPositionIndices(positions, {i, j, k, l})
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
  const std::array<AtomIndex, 4>& indices
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
  const AtomIndex i,
  const AtomIndex j,
  const AtomIndex k,
  const AtomIndex l
) {
  assert(
    validPositionIndices(positions, {i, j, k, l})
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
  ) / 6.0;
}

double getSignedVolume(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndex, 4>& indices
) {
  return getSignedVolume(
    positions,
    indices[0],
    indices[1],
    indices[2],
    indices[3]
  );
}

Eigen::Vector3d averagePosition(
  const Delib::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
) {
  assert(
    validPositionIndices(positions, indices)
  );

  if(indices.size() == 1) {
    return positions[indices.front()].asEigenVector();
  }

  Eigen::Vector3d mean;
  mean.setZero();

  for(const auto& index : indices) {
    mean += positions[index].asEigenVector();
  }

  mean /= indices.size();

  return mean;
}

double distance(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j
) {
  return (i-j).norm();
}

double angle(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k
) {
  Eigen::Vector3d a = i - j;
  Eigen::Vector3d b = k - j;

  return std::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

double dihedral(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
) {
  Eigen::Vector3d a = j - i;
  Eigen::Vector3d b = k - j;
  Eigen::Vector3d c = l - k;

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

double adjustedSignedVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
) {
  return (i - l).dot(
    (j - l).cross(k - l)
  );
}

} // namespace DelibHelpers

} // namespace molassembler
