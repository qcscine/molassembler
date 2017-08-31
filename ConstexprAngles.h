#ifndef INCLUDE_SYMMETRIES_CONSTEXPR_ANGLES_H
#define INCLUDE_SYMMETRIES_CONSTEXPR_ANGLES_H

#include "constexpr_magic/UpperTriangularMatrix.h"

/*! @file
 *
 * Facilitates the calculation of all possible angles between ligand positions
 * of symmetries with many binding sites in a constexpr manner from the
 * coordinates.
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
  return { makeElement(positions, Inds)... };
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

// Binding site positions of square antiprismatic symmetry
constexpr std::array<
  ConstexprMagic::Vector,
  8
> squareAntiprismaticPositions {{
  {{-0.00928803, 0.61156848, 0.79113698}},
  {{0.79562737, 0.60564101, -0.01326839}},
  {{0.79562737, -0.60564101, -0.01326839}},
  {{-0.00928803, -0.61156848, 0.79113698}},
  {{-0.3961716, 0.85216935, -0.34184129}},
  {{0.29375817, 0., -0.95587977}},
  {{-0.3961716, -0.85216935, -0.34184129}},
  {{-0.98308669, 0., 0.18314084}}
}};

/*!
 * An upper triangular matrix containing angles between particules i,j in
 * degrees using the square antiprismatic reference coordinates
 */
constexpr auto squareAntiprismaticAngles = ConstexprMagic::makeUpperTriangularMatrix<8>(
  detail::makeArray<8>(squareAntiprismaticPositions)
);

} // namespace Symmetry

#endif
