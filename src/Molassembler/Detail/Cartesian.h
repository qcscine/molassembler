/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief A series of helpers to interface with the Delib library
 */

#ifndef INCLUDE_MOLASSEMBLER_DELIB_HELPERS_H
#define INCLUDE_MOLASSEMBLER_DELIB_HELPERS_H

#include "Utils/Typenames.h"
#include "Molassembler/Types.h"

#include <Eigen/Geometry>

namespace Scine {
namespace Molassembler {
namespace Cartesian {

/* Reimplementation on vector basis alone */
/*! @brief Averages multiple positions
 *
 * @complexity{@math{\Theta(N)}}
 */
Eigen::Vector3d averagePosition(
  const Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
);

/*! @brief Calculates cartesian distance between two positions
 *
 * @complexity{@math{\Theta(1)}}
 */
double distance(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j
);

/*! @brief Calculates the angle in radians between three positions
 *
 * @complexity{@math{\Theta(1)}}
 */
double angle(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k
);

/*! @brief Calculates the dihedral angle between four positions
 *
 * @complexity{@math{\Theta(1)}}
 * \note Resulting dihedrals are distributed on (-M_PI, M_PI].
 */
double dihedral(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
);

/*! @brief Calculates the adjusted signed tetrahedron volume spanned by four positions
 *
 * The adjusted signed tetrahedron volume @math{V'} is defined from the signed
 * tetrahedron volume @math{V} by @math{V' = 6 V}.
 *
 * @complexity{@math{\Theta(1)}}
 */
double adjustedSignedVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
);

//! @brief Finds the plane of best fit
Eigen::Hyperplane<double, 3> planeOfBestFit(
  const Utils::PositionCollection& positions
);

double planeRmsd(
  const Eigen::Hyperplane<double, 3>& plane,
  const Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
);

/*! @brief Fits a plane to indices and calculates its rms deviation
 *
 * @param positions Full set of positions
 * @param indices Positions to plane fit
 *
 * @complexity{Performs a singular value decomposition. At least linear in the
 * number of indices.}
 */
double planeOfBestFitRmsd(
  const Utils::PositionCollection& positions,
  const std::vector<AtomIndex>& indices
);

//! Reduce any angular value to the boundaries [-pi, pi)
double signedDihedralAngle(double radians);

//! Reduce any angular value to the boundaries [0, 2 pi)
double positiveDihedralAngle(double radians);

/*! Absolute difference between dihedral angles in radians
 *
 * Yields values in the [0, pi) range
 */
double dihedralDifference(double a, double b);

/*! Average of dihedral angles in radians
 *
 * If the angles are exactly opposed, returns the smaller angle offset by pi/2.
 *
 * Returns values in the boundaries [-pi, pi).
 */
double dihedralAverage(double a, double b);

} // namespace Cartesian
} // namespace Molassembler
} // namespace Scine
#endif
