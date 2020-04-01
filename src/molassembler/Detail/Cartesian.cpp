/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Detail/Cartesian.h"

#include <Eigen/Geometry>
#include "molassembler/Temple/Functional.h"

#include <array>

namespace Scine {
namespace molassembler {
namespace cartesian {

bool validPositionIndices(
  const Utils::PositionCollection& positions,
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
  const Utils::PositionCollection& positions,
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

double rmsPlaneDeviation(
  const Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
) {
  const unsigned I = indices.size();
  if(I < 4) {
    throw std::runtime_error("Nonsensical to calculate RMS plane deviation for less than four points");
  }

  /* To find the plane of best fit:
   * - subtract the centroid
   * - calculate the singular value decomposition
   * - normal vector is the left singular vector corresponding to the least
   *   singular value
   *
   * We need to transpose so that we can get away with thin U (Eigen requires a
   * dynamic column count for this, which PositionCollection does not fulfill)
   */
  using ThreeByNPositions = Eigen::Matrix<double, 3, Eigen::Dynamic>;
  ThreeByNPositions relevantPositions (3, I);
  for(unsigned i = 0; i < I; ++i) {
    relevantPositions.col(i) = positions.row(indices[i]).transpose();
  }

  // Calculate the centroid and remove it to get an inertial frame
  const Eigen::Vector3d centroid = relevantPositions.rowwise().sum() / relevantPositions.cols();
  relevantPositions = relevantPositions.colwise() - centroid;

  // SVD values are ordered decreasing -> rightmost column of thin U matrix
  Eigen::JacobiSVD<ThreeByNPositions> decomposition {relevantPositions, Eigen::ComputeThinU | Eigen::ComputeThinV};
  const Eigen::Vector3d planeNormal = decomposition.matrixU().rightCols(1);
  const Eigen::Hyperplane<double, 3> plane {planeNormal, Eigen::Vector3d::Zero()};

  // Calculate the RMS
  double sumOfSquares = 0.0;
  for(unsigned i = 0; i < I; ++i) {
    sumOfSquares += std::pow(plane.signedDistance(relevantPositions.col(i)), 2);
  }
  return std::sqrt(sumOfSquares / I);
}

} // namespace cartesian
} // namespace molassembler
} // namespace Scine
