/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Tetrangle distance bounds smoothing
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_TETRANGLE_SMOOTHING_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_TETRANGLE_SMOOTHING_H

#include <Eigen/Core>

namespace Scine {
namespace Molassembler {
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
 * @warning This doesn't work and I don't know why.
 *
 * @return Tetrangle smoothed bounds matrix
 */
unsigned tetrangleSmooth(Eigen::Ref<Eigen::MatrixXd> bounds);

} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine


#endif
