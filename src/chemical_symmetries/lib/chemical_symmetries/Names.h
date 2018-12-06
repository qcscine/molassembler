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

namespace Symmetry {

//! Enumeration of all contained symmetry names
enum class Name : unsigned {
  Linear, // 2
  Bent,
  TrigonalPlanar, // 3
  CutTetrahedral,
  TShaped,
  Tetrahedral, // 4
  SquarePlanar,
  Seesaw,
  TrigonalPyramidal,
  SquarePyramidal, // 5
  TrigonalBiPyramidal,
  PentagonalPlanar,
  Octahedral, // 6
  TrigonalPrismatic,
  PentagonalPyramidal,
  PentagonalBiPyramidal, // 7
  SquareAntiPrismatic // 8
};

//! Total number of contained symmetries
constexpr unsigned nSymmetries = 17;

} // namespace Symmetry

#endif
