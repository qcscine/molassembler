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
  //! See data::Linear
  Linear,
  //! See data::Bent
  Bent,
  //! See data::TrigonalPlanar
  TrigonalPlanar, // 3
  //! See data::CutTetrahedral
  CutTetrahedral,
  //! See data::TShaped
  TShaped,
  //! See data::Tetrahedral
  Tetrahedral, // 4
  //! See data::SquarePlanar
  SquarePlanar,
  //! See data::Seesaw. Also named disphenoidal.
  Seesaw,
  //! See data::TrigonalPyramidal
  TrigonalPyramidal,
  //! See data::SquarePyramidal
  SquarePyramidal, // 5
  //! See data::TrigonalBiPyramidal
  TrigonalBiPyramidal,
  //! See data::PentagonalPlanar
  PentagonalPlanar,
  //! See data::Octahedral
  Octahedral, // 6
  //! See data::TrigonalPrismatic
  TrigonalPrismatic,
  //! See data::PentagonalPyramidal
  PentagonalPyramidal,
  //! See data::PentagonalBiPyramidal
  PentagonalBiPyramidal, // 7
  //! See data::SquareAntiPrismatic
  SquareAntiPrismatic // 8
};

//! Total number of contained symmetries
constexpr unsigned nSymmetries = 17;

} // namespace Symmetry

} // namespace Scine

#endif
