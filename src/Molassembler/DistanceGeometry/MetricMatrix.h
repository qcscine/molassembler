/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Metric matrix class for DG semantics
 *
 * In the Distance Geometry algorithm, a metric matrix is generated from the
 * distance bounds.
 *
 * The metric matrix then offers the functionality to embed itself into four
 * spatial coordinates. The enclosed closely mirror the algorithms described
 * in rough outline in:
 *
 * - Blaney, J. M., & Dixon, J. S. (2007). Distance Geometry in Molecular
 *   Modeling. Reviews in Computational Chemistry, 5, 299–335.
 *   https://doi.org/10.1002/9780470125823.ch6
 *
 * and in more detail in
 *
 * - Crippen, G. M., & Havel, T. F. (1988). Distance geometry and molecular
 *   conformation (Vol. 74). Taunton: Research Studies Press.
 *   (PDF available as download!)
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_METRIC_MATRIX_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_METRIC_MATRIX_H

#include <Eigen/Core>

#include "Molassembler/DistanceGeometry/DistanceGeometry.h"

namespace Scine {

namespace Molassembler {

namespace DistanceGeometry {

class MetricMatrix {
public:
/* Constructors */
  MetricMatrix() = delete;
  explicit MetricMatrix(Eigen::MatrixXd distanceMatrix);

/* Information */
  /*! @brief Nonmodifiable access to underlying matrix
   *
   * @complexity{@math{\Theta(1)}}
   */
  const Eigen::MatrixXd& access() const;

  /*! @brief Embeds metric matrix into four dimensional space
   *
   * Embeds itself into 4D space, returning a dynamically sized Matrix where
   * every column vector is the coordinates of a particle.
   *
   * @note For Molecules of size 20 and lower, employs full diagonalization. If
   * larger, attempts to calculate only the required eigenpairs. If that fails,
   * falls back on full diagonalization.
   */
  Eigen::MatrixXd embed() const;

  /*! @brief Implements embedding employing full diagonalization
   *
   * Uses Eigen's SelfAdjointEigenSolver to fully diagonalize the matrix,
   * calculating all eigenpairs. Then selects the necessary ones from the full
   * set.
   *
   * @complexity{@math{\Theta(9 N^3)} for the eigenvalue decomposition per
   * Eigen's documentation}
   *
   * @note Faster for roughly N < 20
   */
  Eigen::MatrixXd embedWithFullDiagonalization() const;

/* Operators */
  bool operator == (const MetricMatrix& other) const;

private:
/* Underlying matrix representation */
  Eigen::MatrixXd matrix_;

  void constructFromTemporary_(Eigen::MatrixXd&& distances);
};

} // namespace DistanceGeometry

} // namespace Molassembler

} // namespace Scine

#endif
