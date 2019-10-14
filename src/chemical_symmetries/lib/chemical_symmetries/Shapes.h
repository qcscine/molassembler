/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Defines symmetry names and total count
 *
 * Contains only the symmetry names and count for minimal header inclusion
 * necessities in dependencies
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_NAMES_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_NAMES_H

namespace Scine {

/**
 * @brief Symmetry definitions and properties
 */
namespace Symmetry {

/*! @brief Enumeration of all contained symmetry names
 *
 * List is in order of the number of symmetry positions
 */
enum class Shape : unsigned {
  //! See data::Line
  Line,
  //! See data::Bent
  Bent,
  //! See data::EquilateralTriangle
  EquilateralTriangle, // 3
  //! See data::ApicalTrigonalPyramid
  ApicalTrigonalPyramid,
  //! See data::T
  T,
  //! See data::Tetrahedron
  Tetrahedron, // 4
  //! See data::Square
  Square,
  //! See data::Disphenoid
  Disphenoid,
  //! See data::TrigonalPyramid
  TrigonalPyramid,
  //! See data::SquarePyramid
  SquarePyramid, // 5
  //! See data::TrigonalBipyramid
  TrigonalBipyramid,
  //! See data::Pentagon
  Pentagon,
  //! See data::Octahedron
  Octahedron, // 6
  //! See data::TrigonalPrism
  TrigonalPrism,
  //! See data::PentagonalPyramid
  PentagonalPyramid,
  //! See data::PentagonalBipyramid
  PentagonalBipyramid, // 7
  //! See data::SquareAntiprism
  SquareAntiprism // 8
};

//! Total number of shapes
constexpr unsigned nShapes = 17;
static_assert(nShapes == static_cast<unsigned>(Shape::SquareAntiprism) + 1, "Miscounted?");

} // namespace Symmetry

} // namespace Scine

#endif
