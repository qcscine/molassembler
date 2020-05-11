/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Inertial moments analysis and coordinate reorientation
 */

#include <Eigen/Core>
#include "molassembler/Export.h"

namespace Scine {
namespace Shapes {

using InertialPositionsType = Eigen::Matrix<double, 3, Eigen::Dynamic>;

//! Inertial moments data struct
struct MASM_EXPORT InertialMoments {
  //! Moments values
  Eigen::Vector3d moments;
  //! Moment axes (column-wise)
  Eigen::Matrix3d axes;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/*! @brief Determine the inertial moments of a set of positions
 *
 * @pre Assumes @p positions are in an inertial frame (COM is origin)
 */
MASM_EXPORT InertialMoments principalInertialMoments(const InertialPositionsType& positions);

/**
 * @brief What kind of top is the particle collection?
 */
enum class MASM_EXPORT Top {
  //! Line top: 0 ≅ IA << IB = IC
  Line,
  //! Asymmetric top: IA < IB < IC, degeneracy 0
  Asymmetric,
  //! Prolate top (think rugby football): IA < IB = IC
  Prolate,
  //! Oblate top (think disc): IA = IB < IC
  Oblate,
  //! Spherical top: IA = IB = IC
  Spherical
};

/**
 * @brief Identifies the top of a set of positions and reorients the particle
 *   positions, aligning the main axis along z
 *
 * @param normalizedPositions Particle positions
 *
 * @return The top of the molecule
 */
MASM_EXPORT Top standardizeTop(Eigen::Ref<InertialPositionsType> normalizedPositions);

/**
 * @brief Searches for Cn axes along the coordinate system axes, aligns the
 *   highest order Cn axis found along the z axis
 *
 * @param normalizedPositions Particle positions
 *
 * @note Call standardizeTop() and use this if you get an asymmetric top.
 *
 * @return The order of the highest Cn axis found. If no axis is found along
 *   the coordinate system axes, the result is 1, and particle positions are
 *   unaffected.
 */
MASM_EXPORT unsigned reorientAsymmetricTop(Eigen::Ref<InertialPositionsType> normalizedPositions);

} // namespace Shapes
} // namespace Scine
