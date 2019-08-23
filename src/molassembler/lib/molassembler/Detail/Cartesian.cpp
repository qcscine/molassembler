/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Detail/Cartesian.h"

#include <Eigen/Geometry>
#include "temple/Functional.h"

#include <array>

namespace Scine {

namespace molassembler {

namespace cartesian {

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

} // namespace cartesian

} // namespace molassembler

} // namespace Scine
