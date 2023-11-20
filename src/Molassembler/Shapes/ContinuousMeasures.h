/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Continuous symmetry and shape measures
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_CONTINUOUS_MEASURES_H
#define INCLUDE_MOLASSEMBLER_SHAPES_CONTINUOUS_MEASURES_H

#include "Molassembler/Shapes/PointGroupElements.h"
#include "Molassembler/Shapes/Data.h"

namespace Scine {
namespace Molassembler {
namespace Shapes {

//! @brief Symmetry element, point group, and polyhedral shape continuous metrics
namespace Continuous {

using PositionCollection = Eigen::Matrix<double, 3, Eigen::Dynamic>;

/*! @brief Normalize positions for continuous symmetry measure analysis
 *
 * Reframes to center of mass frame (although no masses exist) and rescales
 * vectors so that the maximum distance is 1)
 */
MASM_EXPORT PositionCollection normalize(const PositionCollection& positions);

//! @brief Continuous symmetry measures for fixed symmetry elements
namespace Fixed {

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
MASM_EXPORT double element(
  const PositionCollection& normalizedPositions,
  const Elements::Rotation& rotation
);

/**
 * @brief Returns the CSM for a fixed reflection symmetry element
 *
 * @param normalizedPositions Particle positions
 * @param reflection Symmetry element of reflection
 *
 * @return The CSM along the reflection plane
 */
MASM_EXPORT double element(
  const PositionCollection& normalizedPositions,
  const Elements::Reflection& reflection
);

/**
 * @brief Returns the CSM for a fixed-axis infinite order rotation axis
 *
 * @param normalizedPositions Particle positions
 * @param axis Axis of the infinite order rotation
 *
 * @return The CSM along the fixed axis
 */
MASM_EXPORT double Cinf(
  const PositionCollection& normalizedPositions,
  const Eigen::Vector3d& axis
);

} // namespace Fixed

/**
 * @brief Optimizes the axis of a rotational symmetry element and calculates the
 *   continuous symmetry measure
 */
MASM_EXPORT std::pair<double, Elements::Rotation> element(
  const PositionCollection& normalizedPositions,
  Elements::Rotation rotation
);

/**
 * @brief Optimizes the norm of a reflection symmetry element and calculates the
 *   continuous symmetry measure
 */
MASM_EXPORT std::pair<double, Elements::Reflection> element(
  const PositionCollection& normalizedPositions,
  Elements::Reflection reflection
);

/*! @brief Calculates the CSM for centroid inversion
 *
 * @note An inversion element cannot be optimized.
 */
MASM_EXPORT double element(
  const PositionCollection& normalizedPositions,
  const Elements::Inversion& /* inversion */
);

/*! @brief Calculates the continuous symmetry measure for an infinite order rotation axis
 */
MASM_EXPORT double Cinf(const PositionCollection& normalizedPositions);

/** @brief Calculates the continuous symmetry measure for a set of particles
 *   and a particular point group
 *
 * @param normalizedPositions
 * @param pointGroup
 *
 * @note This function isn't super stable. The Nelder-Mead over SO(3) works
 * okay, but the initial simplex is a problem and can lead the algorithm into
 * a local minimum. Not very performant either, it would be better with
 * quaternions instead of rotation matrices.
 *
 * @return The continuous symmetry measure
 */
MASM_EXPORT double pointGroup(
  const PositionCollection& normalizedPositions,
  PointGroup pointGroup
);

//! Result of a continuous shape measure calculation
struct MASM_EXPORT ShapeResult {
  //! Lowest value mapping from position indices to shape indices
  std::vector<Vertex> mapping;
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
MASM_EXPORT ShapeResult shapeFaithfulPaperImplementation(
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
MASM_EXPORT ShapeResult shapeAlternateImplementation(
  const PositionCollection& normalizedPositions,
  Shape shape
);

//! Like shapeAlternateImplementation, but the centroid is the last position
MASM_EXPORT ShapeResult shapeAlternateImplementationCentroidLast(
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
MASM_EXPORT ShapeResult shapeHeuristics(
  const PositionCollection& normalizedPositions,
  Shape shape
);

//! @brief Same as shapeHeuristics(), except with set centroid mapping, so faster
MASM_EXPORT ShapeResult shapeHeuristicsCentroidLast(
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
MASM_EXPORT ShapeResult shape(
  const PositionCollection& normalizedPositions,
  Shape shape
);

//! @brief Same as shape(), except with set centroid mapping
MASM_EXPORT ShapeResult shapeCentroidLast(
  const PositionCollection& normalizedPositions,
  Shape shape
);

/*! @brief Calculates minimum distortion angle in radians for shapes A and B
 *
 * Calculates @math{\theta_AB} in:
 *
 * @math{k_XY = \sqrt{\textrm{CShM}A_(B)} = \sqrt{\textrm{CShM}B_(A)} = 10 \sin(\theta_AB)}
 *
 * @warning This function calls shape(), where heuristics are used for
 * particular shape sizes.
 *
 * @complexity{One continuous shape calculation.}
 */
MASM_EXPORT double minimumDistortionAngle(Shape a, Shape b);

/*! @brief Calculates deviation of positions from minimal distortion path between two shapes
 *
 * Calculates
 * @math{\Delta_AB = \frac{1}{\theta_AB}\left[
 *   \arcsin\frac{\sqrt{\textrm{CShM}A_(X)}}{10}
 *   + \arcsin\frac{\sqrt{\textrm{CShM}B_(X)}}{10}
 * \right] - 1} where @math{\theta_AB} is the minimum distortion angle for the
 * shape pair and @math{\textrm{CShM}i_(x)} is the continuous shape measure of
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
MASM_EXPORT double minimalDistortionPathDeviation(
  const PositionCollection& positions,
  Shape a,
  Shape b,
  double minimumDistortionAngle
);

/*! @overload
 *
 * @complexity{Three continuous shape calculations.}
 */
MASM_EXPORT double minimalDistortionPathDeviation(const PositionCollection& positions, Shape a, Shape b);

/*! @brief Beta distribution parameters of shape measures of random point clouds
 *
 * The random point clouds are generated by a zero vector representing the
 * centroid, and points of uniformly distributed direction and normally
 * distributed (mu = 1, stddev = 0.2) length.
 *
 * @param shape Shape for which to determine distribution parameters
 * @param N Number of point clouds to generate
 * @param seed Seeding of the PRNG used for randomness
 *
 * @complexity{N continuous shape calculations, which are basically factorial
 * in the size of the shape}
 *
 * @returns an array of parameters, signifying alpha, beta, loc = 0 and scale
 *   parameters
 */
MASM_EXPORT std::array<double, 4> randomCloudDistributionParameters(
  Shape shape,
  unsigned N,
  unsigned seed
);

/*! @brief Probability that shape measure in set of measures for random point clouds
 *
 * The random point clouds are generated by a zero vector representing the
 * centroid, and points of uniformly distributed direction and normally
 * distributed (mu = 1, stddev = 0.2) length.
 *
 * A beta distribution is fitted against two hundred such continuous shape
 * measures with respect to the chosen shape.
 *
 * The probability is then the cumulative distribution function value of the
 * passed measure in the fitted beta distribution.
 *
 * The distribution parameters for N = 200 random point clouds are precomputed
 * for shapes of sizes <= 8 and hardcoded into the binary. If no beta
 * distribution parameters are available for a shape, returns None.
 *
 * @complexity{Constant}
 */
MASM_EXPORT boost::optional<double> probabilityRandomCloud(double measure, Shape shape);

} // namespace Continuous
} // namespace Shapes
} // namespace Molassembler
} // namespace Scine

#endif
