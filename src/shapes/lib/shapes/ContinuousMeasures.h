/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Continuous symmetry and shape measures
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_CONTINUOUS_MEASURES_H
#define INCLUDE_MOLASSEMBLER_SHAPES_CONTINUOUS_MEASURES_H

#include "shapes/PointGroupElements.h"
#include "shapes/Shapes.h"

namespace Scine {
namespace shapes {

//! @brief Symmetry element, point group, and polyhedral shape continuous metrics
namespace continuous {

using PositionCollection = Eigen::Matrix<double, 3, Eigen::Dynamic>;

/*! @brief Normalize positions for continuous symmetry measure analysis
 *
 * Reframes to center of mass frame (although no masses exist) and rescales
 * vectors so that the maximum distance is 1)
 */
PositionCollection normalize(const PositionCollection& positions);

//! @brief Continuous symmetry measures for fixed symmetry elements
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

/**
 * @brief Optimizes the axis of a rotational symmetry element and calculates the
 *   continuous symmetry measure
 */
std::pair<double, elements::Rotation> element(
  const PositionCollection& normalizedPositions,
  elements::Rotation rotation
);

/**
 * @brief Optimizes the norm of a reflection symmetry element and calculates the
 *   continuous symmetry measure
 */
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
  PointGroup pointGroup
);

//! Result of a continuous shape measure calculation
struct ShapeResult {
  //! Lowest value mapping from position indices to shape indices
  std::vector<unsigned> mapping;
  //! Continuous shape measure value
  double measure;
};

/**
 * @brief Faithful implementation of the continuous shape measure calculation
 *   algorithm from the paper
 *
 * @param normalizedPositions
 * @param shape
 *
 * @complexity{@math{\Theta(N!)}. Note that @math{N} is the size of the shape
 * plus one since a centroid is involved as well.}
 *
 * @return
 */
ShapeResult shapeFaithfulPaperImplementation(
  const PositionCollection& normalizedPositions,
  Shape shape
);

/**
 * @brief Slightly optimized implmentation of the continuous shape measure
 *   calculation algorithm
 *
 * Paper says:
 * - For each index mapping
 *   - Minimize over rotation
 *   - Minimize over isotropic scaling factor.
 *
 * This does:
 * - For each index mapping
 *   - Minimize over rotation
 * - For the best index mapping, minimize over isotropic scaling factor
 *
 * and is hence faster. Gives the same results.
 *
 * @param normalizedPositions
 * @param shape
 *
 * @return
 */
ShapeResult shapeAlternateImplementation(
  const PositionCollection& normalizedPositions,
  Shape shape
);

//! Like shapeAlternateImplementation, but the centroid is the last position
ShapeResult shapeAlternateImplementationCentroidLast(
  const PositionCollection& normalizedPositions,
  Shape shape
);

/**
 * @brief Calculates the continuous shape measure of a set of coordinates with
 *   respect to a particular shape using heuristics
 *
 * Since shapeFaithfulPaperImplementation() scales as @math{\Theta(N!)} and we
 * have shapes with 12 vertices plus the centroid, we need heuristics to speed
 * up the shape calculation while not sacrificing too much accuracy.
 *
 * Heuristic used here: For all tuples of five positions, align the positions,
 * then greedily choose the best next sequence alignments until all positions
 * are matched.
 *
 * Judgement of accuracy / speed tradeoff is that this is worth using from size
 * 6 onwards in debug builds and from size 9 onwards in release builds.
 *
 * @param normalizedPositions set of coordinates to compare with the shape
 * @param shape Reference shape to compare against
 *
 * @complexity{Approximately @math{\Omega(\frac{N!}{(N-5)!})} quaternion fits,
 * where @math{N} is the number of positions being matched. Note that @math{N}
 * is typically the size of the shape plus one since a centroid is involved as
 * well. For N > 12, other terms may dominate complexity, but this is yet
 * untested.}
 *
 * @note Works well for positions deviating little from the ideal shape. For
 * large deviations, exhibits small relative errors. To explore the
 * characteristics of this function, you can try out the shape analysis binary.
 *
 * @return The continuous shape measure determined by heuristics
 */
ShapeResult shapeHeuristics(
  const PositionCollection& normalizedPositions,
  Shape shape
);

//! @brief Same as shapeHeuristics(), except with set centroid mapping, so faster
ShapeResult shapeHeuristicsCentroidLast(
  const PositionCollection& normalizedPositions,
  Shape shape
);

/**
 * @brief Forwarding function to calculate the continuous shape measure
 *
 * Forwards its call to shapeAlternateImplementation() by default. In debug
 * builds, forwards its call to shapeHeuristics() from shape size 6 onwards. In
 * release builds, forwards its call to shapeHeuristics() from shape size 9
 * onwards.
 */
ShapeResult shape(
  const PositionCollection& normalizedPositions,
  Shape shape
);

//! @brief Same as shape(), except with set centroid mapping
ShapeResult shapeCentroidLast(
  const PositionCollection& normalizedPositions,
  const Shape shape
);

/*! @brief Calculates minimum distortion angle in radians for shapes A and B
 *
 * Calculates @math{\theta_AB} in:
 *
 * @math{k_XY = \sqrt{\textrm{CShM}_A(B)} = \sqrt{\textrm{CShM}_B(A)} = 10 \sin(\theta_AB)}
 *
 * @warning This function calls shape(), where heuristics are used for
 * particular shape sizes.
 *
 * @complexity{One continuous shape calculation.}
 */
double minimumDistortionAngle(Shape a, Shape b);

/*! @brief Calculates deviation of positions from minimal distortion path between two shapes
 *
 * Calculates
 * @math{\Delta_AB = \frac{1}{\theta_AB}\left[
 *   \arcsin\frac{\sqrt{\textrm{CShM}_A(X)}}{10}
 *   + \arcsin\frac{\sqrt{\textrm{CShM}_B(X)}}{10}
 * \right] - 1} where @math{\theta_AB} is the minimum distortion angle for the
 * shape pair and @math{\textrm{CShM}_i(x)} is the continuous shape measure of
 * the positions @math{x} with regards to the shape @math{i}.
 *
 * This function form avoids two shape calculations if the minimum distortion
 * angle between a and b is known.
 *
 * @warning This function calls shape(), where heuristics are used for
 * particular shape sizes.
 *
 * @complexity{Two continuous shape calculations.}
 */
double minimalDistortionPathDeviation(
  const PositionCollection& positions,
  Shape a,
  Shape b,
  const double minimumDistortionAngle
);

/*! @overload
 *
 * @complexity{Three continuous shape calculations.}
 */
double minimalDistortionPathDeviation(const PositionCollection& positions, Shape a, Shape b);


} // namespace continuous
} // namespace shapes
} // namespace Scine

#endif
