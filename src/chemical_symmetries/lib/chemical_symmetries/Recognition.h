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

struct InertialMoments {
  Eigen::Vector3d moments;
  Eigen::Matrix3d axes;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

//! @pre Assumes this is an inertial frame (COM is origin)
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

Top standardizeTop(Eigen::Ref<PositionCollection> normalizedPositions);

PointGroup flowchart(
  const PositionCollection& normalizedPositions,
  const Top top
);

/**
 * @brief Namespace for calculation of continuous symmetry measures
 */
namespace csm {

/*! @brief Minimizes CSM for a point group, case: G = P
 *
 * This minimizes the continuous symmetry measure for the case that the number
 * of group symmetry elements matches the number of particles.
 */
double allSymmetryElements(
  const PositionCollection& normalizedPositions,
  const std::vector<std::unique_ptr<elements::SymmetryElement>>& elements,
  std::vector<unsigned> particleIndices
);

/*! @brief Minimizes CSM for a point group, case G = l * P
 *
 * This minimizes the continuous symmetry measure for the case that the number
 * of group symmetry elements is a multiple l of the number of particles.
 *
 * @todo consider particleIndices by-ref
 */
double groupedSymmetryElements(
  const PositionCollection& normalizedPositions,
  std::vector<unsigned> particleIndices,
  const std::vector<std::unique_ptr<elements::SymmetryElement>>& elements,
  const elements::ElementGrouping& elementGrouping
);

/** @brief Calculates the continuous symmetry measure for a set of particles
 *   and a particular point group
 *
 * @param normalizedPositions
 * @param pointGroup
 *
 * @return The continuous symmetry measure
 */
double point_group(
  const PositionCollection& normalizedPositions,
  const PointGroup pointGroup
);

} // namespace csm

namespace detail {

/*! @brief Normalize positions for continuous symmetry measure analysis
 *
 * Reframes to center of mass frame (although no masses exist) and rescales
 * vectors so that the maximum distance is 1)
 */
PositionCollection normalize(const PositionCollection& positions);

} // namespace detail

} // namespace Symmetry
} // namespace Scine

#endif
