/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Continuous symmetry and shape measures
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_CONTINUOUS_MEASURES_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_CONTINUOUS_MEASURES_H

#include "chemical_symmetries/PointGroupElements.h"
#include "chemical_symmetries/Shapes.h"

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

/* Regarding CShapeM minimization, there is a bit of a conundrum: The paper says
 * - For each index mapping
 *   - Minimize over rotation
 *   - Minimize over scaling
 *
 * But it's a lot faster to
 * - Minimize over rotation (while minimizing CShapeM over all permutations)
 * - Minimize over scaling (using best permutation from rotation step)
 *
 * And to me it's not immediately apparent why this should be worse, especially
 * considering that minimizing over permutations while calculating CShapeM
 * should be smooth. Maybe it has local minima that the paper procedure wouldn't?
 *
 * But it's odd. The paper suggests pre-pairing off vertices to reduce cost "in
 * most cases". Not sure which is more dangerous (pre-pairing prior to any
 * minimizations or reversing minimization order and reusing pairing).
 *
 * Furthermore, both variants get different results despite fitting looking
 * veeery similar. Messing with minimization epsilons doesn't get me anywhere.
 *
 * Another way to do shape minimization without full permutational work might
 * be by sequential alignment: Manage an unordered map keeping track of the
 * emerging index permutation. Always add that mapping with the minimal cost to
 * the map and realign by orientation minimization.
 */

/**
 * @brief Calculates the continuous shape measure of a set of coordinates with
 *   respect to a particular shape
 *
 * @param normalizedPositions
 * @param shape
 *
 * @throws std::logic_error If the number of passed positions does not match
 * the size of the shape
 *
 * @return The continuous shape measure
 */
double shape(
  const PositionCollection& normalizedPositions,
  const Shape shape
);

double shapeFaithfulPaperImplementation(
  const PositionCollection& normalizedPositions,
  const Shape shape
);

double shapeAlternateImplementation(
  const PositionCollection& normalizedPositions,
  const Shape shape
);

} // namespace continuous
} // namespace Symmetry
} // namespace Scine

#endif
