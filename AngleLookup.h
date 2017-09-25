#ifndef INCLUDE_SYMMETRIES_CONSTEXPR_ANGLES_H
#define INCLUDE_SYMMETRIES_CONSTEXPR_ANGLES_H

#include "constexpr_magic/UpperTriangularMatrix.h"

/*! @file
 *
 * Facilitates the calculation of all possible angles between ligand positions
 * of symmetries with many binding sites in a constexpr manner from the
 * coordinates.
 */

/*! @file 
 *
 * Provides functionality for converting a set of 3D atom positions into an
 * array of all possible angles between the positions in a \c constexpr fashion.
 *
 * Explicitly creates this array for a set of atoms arranged in a square
 * antiprismatic geometry.
 */

namespace Symmetry {

namespace detail {

template<unsigned long size>
constexpr double makeElement(
  const std::array<ConstexprMagic::Vector, size>& positions,
  const size_t& i
) {
  // Get i-j matrix indices from the linear index
  const auto indexPair = ConstexprMagic::UpperTriangularMatrixImpl::index_conversion::toDoubleIndex<size>(i);

  // Calculate the angle
  return ConstexprMagic::angle(
    positions[indexPair.first],
    positions[indexPair.second]
  );
}

template<unsigned long size, size_t... Inds>
constexpr std::array<double, size * (size - 1) / 2> makeArrayImpl(
  const std::array<ConstexprMagic::Vector, size>& positions,
  std::integer_sequence<size_t, Inds...>
) {
  // Expand the parameter pack for each individual linear index
  return {{ makeElement(positions, Inds)... }};
}

/* Entry point for array creation, calculates the required dimension of the
 * upper triangular matrix required to store all angles
 */
template<unsigned long size>
constexpr std::array<double, size * (size - 1) / 2> makeArray(
  const std::array<ConstexprMagic::Vector, size>& positions
) {
  return makeArrayImpl(
    positions,
    std::make_index_sequence<size * (size - 1) / 2>{}
  );
}

} // namespace detail

} // namespace Symmetry

#endif
