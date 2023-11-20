/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Detail/Cartesian.h"

#include "Molassembler/Temple/Functional.h"

#include <array>

namespace Scine {
namespace Molassembler {
namespace Cartesian {

bool validPositionIndices(
  const Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
) {
  const unsigned nRows = positions.rows();
  return Temple::all_of(
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

Eigen::Hyperplane<double, 3> planeOfBestFit(const Utils::PositionCollection& positions) {
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
  ThreeByNPositions transpose = positions.transpose();

  // Calculate the centroid and remove it to get an inertial frame
  const Eigen::Vector3d centroid = transpose.rowwise().sum() / transpose.cols();
  transpose = transpose.colwise() - centroid;

  // SVD values are ordered decreasing -> rightmost column of thin U matrix
  Eigen::JacobiSVD<ThreeByNPositions> decomposition {transpose, Eigen::ComputeThinU | Eigen::ComputeThinV};
  return {decomposition.matrixU().rightCols(1), centroid};
}

double planeRmsd(
  const Eigen::Hyperplane<double, 3>& plane,
  const Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
) {
  const unsigned I = indices.size();
  double sumOfSquares = 0.0;
  for(unsigned i = 0; i < I; ++i) {
    sumOfSquares += std::pow(plane.signedDistance(positions.row(indices[i])), 2);
  }
  return std::sqrt(sumOfSquares / I);
}

double planeOfBestFitRmsd(
  const Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
) {
  const unsigned I = indices.size();
  if(I < 4) {
    throw std::runtime_error("Nonsensical to calculate RMS plane deviation for less than four points");
  }

  // Select the relevant positions only
  Utils::PositionCollection relevantPositions (I, 3);
  for(unsigned i = 0; i < I; ++i) {
    relevantPositions.row(i) = positions.row(indices[i]);
  }

  const Eigen::Hyperplane<double, 3> plane = planeOfBestFit(relevantPositions);
  return planeRmsd(plane, positions, indices);
}

double signedDihedralAngle(const double radians) {
  return radians - std::floor((radians + M_PI) / (2 * M_PI)) * 2 * M_PI;
}

double positiveDihedralAngle(const double radians) {
  return radians - std::floor(radians / (2 * M_PI)) * 2 * M_PI;
}

double dihedralDifference(const double a, const double b) {
  double midpointComparison = std::fabs(signedDihedralAngle(a) - signedDihedralAngle(b));
  double zeroboundedComparison = std::fabs(positiveDihedralAngle(a) - positiveDihedralAngle(b));
  return std::min(midpointComparison, zeroboundedComparison);
}

double dihedralAverage(const double a, const double b) {
  const double s = (std::sin(a) + std::sin(b)) / 2;
  const double c = (std::cos(a) + std::cos(b)) / 2;

  // What to do when angles are opposed (offset by pi)
  if(s * s + c * c <= 1e-20) {
    return signedDihedralAngle(std::min(a, b) + M_PI / 2);
  }

  return std::atan2(s, c);
}

} // namespace Cartesian
} // namespace Molassembler
} // namespace Scine
