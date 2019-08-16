/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Detail/DelibHelpers.h"

#include <Eigen/Geometry>
#include "temple/Functional.h"

#include <array>

namespace Scine {

namespace molassembler {

namespace DelibHelpers {

bool validPositionIndices(
  const Scine::Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
) {
  const unsigned nRows = positions.rows();
  return temple::all_of(
    indices,
    [&](const AtomIndex i) -> bool {
      return i < nRows;
    }
  );
}

double getDihedral(
  const Scine::Utils::PositionCollection& positions,
  const AtomIndex i,
  const AtomIndex j,
  const AtomIndex k,
  const AtomIndex l
) {
  assert(
    validPositionIndices(positions, {i, j, k, l})
  );

  Eigen::Vector3d a = positions.row(j) - positions.row(i);
  Eigen::Vector3d b = positions.row(k) - positions.row(j);
  Eigen::Vector3d c = positions.row(l) - positions.row(k);

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
  const Scine::Utils::PositionCollection& positions,
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
  const Scine::Utils::PositionCollection& positions,
  const AtomIndex i,
  const AtomIndex j,
  const AtomIndex k,
  const AtomIndex l
) {
  assert(
    validPositionIndices(positions, {i, j, k, l})
  );

  return (
    positions.row(i)
    - positions.row(l)
  ).dot(
    (
      positions.row(j)
      - positions.row(l)
    ).cross(
      positions.row(k)
      - positions.row(l)
    )
  ) / 6.0;
}

double getSignedVolume(
  const Scine::Utils::PositionCollection& positions,
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
  const Scine::Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
) {
  assert(!indices.empty());
  assert(
    validPositionIndices(positions, indices)
  );

  if(indices.size() == 1) {
    return positions.row(indices.front());
  }

  Eigen::Vector3d mean;
  mean.setZero();

  for(const auto& index : indices) {
    mean += positions.row(index);
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

} // namespace Scine
