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
  /* 2 */
  //! See data::Line
  Line,
  //! See data::Bent
  Bent,

  /* 3 */
  //! See data::EquilateralTriangle
  EquilateralTriangle,
  //! See data::VacantTetrahedron
  VacantTetrahedron,
  //! See data::T
  T,

  /* 4 */
  //! See data::Tetrahedron
  Tetrahedron,
  //! See data::Square
  Square,
  //! See data::Seesaw
  Seesaw,
  //! See data::TrigonalPyramid
  TrigonalPyramid,

  /* 5 */
  //! See data::SquarePyramid
  SquarePyramid,
  //! See data::TrigonalBipyramid
  TrigonalBipyramid,
  //! See data::Pentagon
  Pentagon,

  /* 6 */
  //! See data::Octahedron
  Octahedron, // 6
  //! See data::TrigonalPrism
  TrigonalPrism,
  //! See data::PentagonalPyramid
  PentagonalPyramid,
  //! See data::Hexagon
  Hexagon,

  /* 7 */
  //! See data::PentagonalBipyramid
  PentagonalBipyramid,
  //! See data::CappedOctahedron
  CappedOctahedron,
  //! See data::CappedTrigonalPrism
  CappedTrigonalPrism,

  /* 8 */
  //! See data::SquareAntiprism
  SquareAntiprism,
  //! See data::Cube
  Cube,
  //! See data::TrigonalDodecahedron
  TrigonalDodecahedron,
  //! See data::HexagonalBipyramid
  HexagonalBipyramid,

  /* 9 */
  //! See data::TricappedTrigonalPrism
  TricappedTrigonalPrism,
  //! See data::CappedSquareAntiPrism
  CappedSquareAntiprism,
  //! See data::HeptagonalBipyramid
  HeptagonalBipyramid,

  /* 10 */
  //! See data::BicappedSquareAntiprism
  BicappedSquareAntiprism,

  /* 11 */
  //! See data::EdgeContractedIcosahedron
  EdgeContractedIcosahedron,

  /* 12 */
  //! See data::Icosahedron
  Icosahedron,
  //! See data::Cuboctahedron
  Cuboctahedron
};

//! Total number of shapes
constexpr unsigned nShapes = 30;
static_assert(nShapes == static_cast<unsigned>(Shape::Cuboctahedron) + 1, "Miscounted?");

} // namespace Symmetry

} // namespace Scine

#endif
