/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Tetrangle distance bounds smoothing
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_TETRANGLE_SMOOTHING_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_TETRANGLE_SMOOTHING_H

#include <Eigen/Core>

namespace Scine {
namespace molassembler {
namespace DistanceGeometry {

/**
 * @brief Smoothes the bounds matrix using tetrangle inequalities
 *
 * @pre Assumes @p bounds has already been triangle inequality smoothed
 *
 * @param bounds A square matrix with zeros on the diagonal and
 *   B(j, i) <= B(i, j) for all i < j (lower bounds on strictly lower triangle,
 *   upper bounds on strictly upper triangle)
 *
 * @complexity{@math{\Theta(N^4)}}
 *
 * @return Tetrangle smoothed bounds matrix
 */
Eigen::MatrixXd tetrangleSmooth(Eigen::MatrixXd bounds);

} // namespace DistanceGeometry
} // namespace molassembler
} // namespace Scine


#endif
