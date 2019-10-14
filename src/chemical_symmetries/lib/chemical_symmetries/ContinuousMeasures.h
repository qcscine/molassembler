/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Continuous symmetry and shape measures
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_CONTINUOUS_MEASURES_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_CONTINUOUS_MEASURES_H

#include "chemical_symmetries/PointGroupElements.h"

namespace Scine {
namespace Symmetry {
namespace continuous {

using PositionCollection = Eigen::Matrix<double, 3, Eigen::Dynamic>;

/*! @brief Normalize positions for continuous symmetry measure analysis
 *
 * Reframes to center of mass frame (although no masses exist) and rescales
 * vectors so that the maximum distance is 1)
 */
PositionCollection normalize(const PositionCollection& positions);

namespace fixed {

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

/**
 * @brief Returns the CSM for a fixed reflection symmetry element
 *
 * @param normalizedPositions Particle positions
 * @param reflection Symmetry element of reflection
 *
 * @return The CSM along the reflection plane
 */
double element(
  const PositionCollection& normalizedPositions,
  const elements::Reflection& reflection
);

/**
 * @brief Returns the CSM for a fixed-axis infinite order rotation axis
 *
 * @param normalizedPositions Particle positions
 * @param axis Axis of the infinite order rotation
 *
 * @return The CSM along the fixed axis
 */
double Cinf(
  const PositionCollection& normalizedPositions,
  const Eigen::Vector3d& axis
);

} // namespace fixed

std::pair<double, elements::Rotation> element(
  const PositionCollection& normalizedPositions,
  elements::Rotation rotation
);

std::pair<double, elements::Reflection> element(
  const PositionCollection& normalizedPositions,
  elements::Reflection reflection
);

/*! @brief Calculates the CSM for centroid inversion
 *
 * @note An inversion element cannot be optimized.
 */
double element(
  const PositionCollection& normalizedPositions,
  const elements::Inversion& /* inversion */
);

/*! @brief Calculates the continuous symmetry measure for an infinite order rotation axis
 */
double Cinf(const PositionCollection& normalizedPositions);

/** @brief Calculates the continuous symmetry measure for a set of particles
 *   and a particular point group
 *
 * @param normalizedPositions
 * @param pointGroup
 *
 * @return The continuous symmetry measure
 */
double pointGroup(
  const PositionCollection& normalizedPositions,
  const PointGroup pointGroup
);

} // namespace continuous
} // namespace Symmetry
} // namespace Scine

#endif
