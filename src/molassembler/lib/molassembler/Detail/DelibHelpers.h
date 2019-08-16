/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief A series of helpers to interface with the Delib library
 */

#ifndef INCLUDE_MOLASSEMBLER_DELIB_HELPERS_H
#define INCLUDE_MOLASSEMBLER_DELIB_HELPERS_H

#include "Utils/Typenames.h"
#include "molassembler/Types.h"

namespace Scine {

namespace molassembler {

namespace DelibHelpers {

/*! @brief Calculates the dihedral angle in radians defined by four positions
 *
 * @complexity{@math{\Theta(1)}}
 * \note Resulting dihedrals are distributed on (-M_PI, M_PI].
 */
double getDihedral(
  const Scine::Utils::PositionCollection& positions,
  AtomIndex i,
  AtomIndex j,
  AtomIndex k,
  AtomIndex l
);

/*! @brief Calculates the dihedral angle in radians defined by four positions
 *
 * @complexity{@math{\Theta(1)}}
 * \note Resulting dihedrals are distributed on (-M_PI, M_PI].
 */
double getDihedral(
  const Scine::Utils::PositionCollection& positions,
  const std::array<AtomIndex, 4>& indices
);

/*! @brief Calculates signed tetrahedron volume spanned by four positions
 *
 * @complexity{@math{\Theta(1)}}
 */
double getSignedVolume(
  const Scine::Utils::PositionCollection& positions,
  AtomIndex i,
  AtomIndex j,
  AtomIndex k,
  AtomIndex l
);

/*! @brief Calculates signed tetrahedron volume spanned by four positions
 *
 * @complexity{@math{\Theta(1)}}
 */
double getSignedVolume(
  const Scine::Utils::PositionCollection& positions,
  const std::array<AtomIndex, 4>& indices
);

/* Reimplementation on vector basis alone */
/*! @brief Averages multiple positions
 *
 * @complexity{@math{\Theta(N)}}
 */
Eigen::Vector3d averagePosition(
  const Scine::Utils::PositionCollection& positions,
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
 */
double adjustedSignedVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
);

} // namespace DelibHelpers

} // namespace molassembler

} // namespace Scine
#endif
