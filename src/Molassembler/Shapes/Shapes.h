/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Defines symmetry names and total count
 *
 * Contains only the symmetry names and count for minimal header inclusion
 * necessities in dependencies
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_NAMES_H
#define INCLUDE_MOLASSEMBLER_SHAPES_NAMES_H

#include "Molassembler/Export.h"

namespace Scine {
namespace Molassembler {

/**
 * @brief Symmetry definitions and properties
 */
namespace Shapes {

/*! @brief Enumeration of all contained symmetry names
 *
 * List is in order of the number of symmetry positions
 */
enum class MASM_EXPORT Shape : unsigned {
  /* 2 */
  //! See Data::Line
  Line,
  //! See Data::Bent
  Bent,

  /* 3 */
  //! See Data::EquilateralTriangle
  EquilateralTriangle,
  //! See Data::VacantTetrahedron
  VacantTetrahedron,
  //! See Data::T
  T,

  /* 4 */
  //! See Data::Tetrahedron
  Tetrahedron,
  //! See Data::Square
  Square,
  //! See Data::Seesaw
  Seesaw,
  //! See Data::TrigonalPyramid
  TrigonalPyramid,

  /* 5 */
  //! See Data::SquarePyramid
  SquarePyramid,
  //! See Data::TrigonalBipyramid
  TrigonalBipyramid,
  //! See Data::Pentagon
  Pentagon,

  /* 6 */
  //! See Data::Octahedron
  Octahedron, // 6
  //! See Data::TrigonalPrism
  TrigonalPrism,
  //! See Data::PentagonalPyramid
  PentagonalPyramid,
  //! See Data::Hexagon
  Hexagon,

  /* 7 */
  //! See Data::PentagonalBipyramid
  PentagonalBipyramid,
  //! See Data::CappedOctahedron
  CappedOctahedron,
  //! See Data::CappedTrigonalPrism
  CappedTrigonalPrism,

  /* 8 */
  //! See Data::SquareAntiprism
  SquareAntiprism,
  //! See Data::Cube
  Cube,
  //! See Data::TrigonalDodecahedron
  TrigonalDodecahedron,
  //! See Data::HexagonalBipyramid
  HexagonalBipyramid,

  /* 9 */
  //! See Data::TricappedTrigonalPrism
  TricappedTrigonalPrism,
  //! See Data::CappedSquareAntiPrism
  CappedSquareAntiprism,
  //! See Data::HeptagonalBipyramid
  HeptagonalBipyramid,

  /* 10 */
  //! See Data::BicappedSquareAntiprism
  BicappedSquareAntiprism,

  /* 11 */
  //! See Data::EdgeContractedIcosahedron
  EdgeContractedIcosahedron,

  /* 12 */
  //! See Data::Icosahedron
  Icosahedron,
  //! See Data::Cuboctahedron
  Cuboctahedron
};

//! Total number of shapes
constexpr unsigned nShapes = 30;
static_assert(nShapes == static_cast<unsigned>(Shape::Cuboctahedron) + 1, "Miscounted?");

} // namespace Shapes
} // namespace Molassembler
} // namespace Scine

#endif
