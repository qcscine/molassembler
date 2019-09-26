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
enum class Name : unsigned {
  /* 2 */
  //! See data::Linear
  Linear,
  //! See data::Bent
  Bent,

  /* 3 */
  //! See data::TrigonalPlanar
  TrigonalPlanar,
  //! See data::CutTetrahedral
  CutTetrahedral,
  //! See data::TShaped
  TShaped,

  /* 4 */
  //! See data::Tetrahedral
  Tetrahedral,
  //! See data::SquarePlanar
  SquarePlanar,
  //! See data::Seesaw.
  Seesaw,
  //! Alternate name for Seesaw
  Disphenoidal = Seesaw,
  //! See data::TrigonalPyramidal
  TrigonalPyramidal,

  /* 5 */
  //! See data::SquarePyramidal
  SquarePyramidal,
  //! See data::TrigonalBiPyramidal
  TrigonalBiPyramidal,
  //! See data::PentagonalPlanar
  PentagonalPlanar,

  /* 6 */
  //! See data::Octahedral
  Octahedral,
  //! See data::TrigonalPrismatic
  TrigonalPrismatic,
  //! See data::PentagonalPyramidal
  PentagonalPyramidal,

  /* 7 */
  //! See data::PentagonalBiPyramidal
  PentagonalBiPyramidal,
  //! See data::CappedOctahedral
  CappedOctahedral,
  //! See data::CappedTrigonalPrismatic
  CappedTrigonalPrismatic,

  /* 8 */
  //! See data::SquareAntiPrismatic
  SquareAntiPrismatic,
  //! See data::Cubic
  Cubic,
  //! See data::BicappedTrigonalPrismatic
  BicappedTrigonalPrismatic,
  //! See data::Dodecahedral
  Dodecahedral,
  //! Alternate name for Dodecahedral
  SnubDisphenoidal = Dodecahedral,
  //! See data::HexagonalBiPyramidal
  HexagonalBiPyramidal,

  /* 9 */
  //! See data::TricappedTrigonalPrismatic
  TricappedTrigonalPrismatic,
  //! See data::CappedSquareAntiPrismatic
  CappedSquareAntiPrismatic,
  //! See data::HeptagonalBiPyramidal
  HeptagonalBiPyramidal,

  /* 10 */
  //! See data::BicappedSquareAntiPrismatic
  BicappedSquareAntiPrismatic,

  /* 12 */
  //! See data::Icosahedral
  Icosahedral,
  //! See data::Cuboctahedron
  Cuboctahedron
};

//! Total number of contained symmetries
constexpr unsigned nSymmetries = 29;
static_assert(nSymmetries == static_cast<unsigned>(Name::Cuboctahedron) + 1, "Miscounted?");

} // namespace Symmetry

} // namespace Scine

#endif
