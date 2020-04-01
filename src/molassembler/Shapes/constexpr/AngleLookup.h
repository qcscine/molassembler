/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Generate lookup table for a symmetry's angles
 *
 * Facilitates the calculation of all possible angles between ligand positions
 * of symmetries with many binding sites in a constexpr manner from the
 * coordinates.
 *
 * Provides functionality for converting a set of 3D atom positions into an
 * array of all possible angles between the positions in a \c constexpr fashion.
 *
 * Explicitly creates this array for a set of atoms arranged in a square
 * antiprismatic geometry.
 */

#ifndef INCLUDE_SYMMETRIES_CONSTEXPR_ANGLES_H
#define INCLUDE_SYMMETRIES_CONSTEXPR_ANGLES_H

#include "molassembler/Temple/constexpr/UpperTriangularMatrix.h"

namespace Scine {

namespace shapes {

namespace detail {

/*! @brief Calculate the i-th element of the angle upper triangular matrix
 *
 * @complexity{@math{\Theta(1)}}
 *
 * @param positions The positions array of the symmetry of interest
 * @param i Single-index index into the linear storage of the upper triangular
 *   matrix
 */
template<unsigned long size>
constexpr double makeElement(
  const std::array<temple::Vector, size>& positions,
  const size_t i
) {
  // Get i-j matrix indices from the linear index
  const auto indexPair = temple::UpperTriangularMatrixImpl::index_conversion::toDoubleIndex<size>(i);

  // Calculate the angle
  return temple::angle(
    positions[indexPair.first],
    positions[indexPair.second]
  );
}

/*! @brief Generate the linear storage of the angle upper triangular matrix
 *
 * @complexity{@math{\Theta(N^2)}}
 *
 * @param positions The positions array of the symmetry of interest
 * @param inds An integer sequence of appropriate length for the desired
 *   symmetry
 */
template<unsigned long size, size_t... Inds>
constexpr std::array<double, size * (size - 1) / 2> makeArrayImpl(
  const std::array<temple::Vector, size>& positions,
  std::integer_sequence<size_t, Inds...> /* inds */
) {
  // Expand the parameter pack for each individual linear index
  return {{ makeElement(positions, Inds)... }};
}

/*! @brief Generate the linear storage of the angle upper triangular matrix
 *
 * @complexity{@math{\Theta(N^2)}}
 *
 * Entry point for array creation, calculates the required dimension of the
 * linear array underlying the upper triangular matrix required to store all
 * angles for a particular symmetry
 *
 * @param positions The positions array of the symmetry of interest
 */
template<unsigned long size>
constexpr std::array<double, size * (size - 1) / 2> makeArray(
  const std::array<temple::Vector, size>& positions
) {
  return makeArrayImpl(
    positions,
    std::make_index_sequence<size * (size - 1) / 2>{}
  );
}

} // namespace detail

} // namespace shapes

} // namespace Scine

#endif
