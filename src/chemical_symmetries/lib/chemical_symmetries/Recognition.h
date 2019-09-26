/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Analyze coordinates for point group symmetry
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_RECOGNITION_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_RECOGNITION_H

#include "chemical_symmetries/PointGroups.h"

namespace Scine {
namespace Symmetry {

using PositionCollection = Eigen::Matrix<double, 3, Eigen::Dynamic>;

//! Inertial moments data struct
struct InertialMoments {
  //! Moments values
  Eigen::Vector3d moments;
  //! Moment axes (column-wise)
  Eigen::Matrix3d axes;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/*! @brief Determine the inertial moments of a set of positions
 *
 * @pre Assumes this is an inertial frame (COM is origin)
 */
InertialMoments principalInertialMoments(
  const PositionCollection& normalizedPositions
);

/**
 * @brief What kind of top is the particle collection?
 */
enum class Top {
  //! Linear top: 0 ≅ IA << IB = IC
  Linear,
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
Top standardizeTop(Eigen::Ref<PositionCollection> normalizedPositions);

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
unsigned reorientAsymmetricTop(Eigen::Ref<PositionCollection> normalizedPositions);

/**
 * @brief Attempts to find the point group of a set of positions by flowcharting
 *   individual symmetry elements.
 *
 * @warning Incomplete
 *
 * @param normalizedPositions Particle positions
 * @param top The result of top standardization
 *
 * @return The point group of the particle positions
 */
PointGroup flowchart(
  const PositionCollection& normalizedPositions,
  Top top
);

/**
 * @brief Namespace for calculation of continuous symmetry measures
 */
namespace csm {

/** @brief Calculates the continuous symmetry measure for a set of particles
 *   and a particular point group
 *
 * @param normalizedPositions
 * @param pointGroup
 *
 * @return The continuous symmetry measure
 */
boost::optional<double> pointGroup(
  const PositionCollection& normalizedPositions,
  const PointGroup pointGroup
);

/**
 * @brief Returns the CSM for a Rotation symmetry element along the rotation
 *   axis without optimizing the coordinates' rotation
 *
 * @param normalizedPositions Particle positions
 * @param rotation Symmetry element of rotation Cn/Sn
 *
 * @pre @p rotation power is one, and its axis is normalized (latter is
 *   guaranteed by its constructor)
 *
 * @return The CSM along the fixed axis of rotation
 */
double element(
  const PositionCollection& normalizedPositions,
  const elements::Rotation& rotation
);

// TODO more elements, optimize elements?

} // namespace csm

namespace detail {

/*! @brief Normalize positions for continuous symmetry measure analysis
 *
 * Reframes to center of mass frame (although no masses exist) and rescales
 * vectors so that the maximum distance is 1)
 */
PositionCollection normalize(const PositionCollection& positions);

PointGroup linear(const PositionCollection& normalizedPositions);

} // namespace detail

} // namespace Symmetry
} // namespace Scine

#endif
