/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Point group enum
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_POINT_GROUPS_H
#define INCLUDE_MOLASSEMBLER_SHAPES_POINT_GROUPS_H

#include "Molassembler/Export.h"

namespace Scine {
namespace Molassembler {
namespace Shapes {

/**
 * @brief Point groups
 */
enum class MASM_EXPORT PointGroup : unsigned {
  C1, Ci, Cs,
  C2, C3, C4, C5, C6, C7, C8,
  C2h, C3h, C4h, C5h, C6h, C7h, C8h,
  C2v, C3v, C4v, C5v, C6v, C7v, C8v,
  S4, S6, S8,
  D2, D3, D4, D5, D6, D7, D8,
  D2h, D3h, D4h, D5h, D6h, D7h, D8h,
  D2d, D3d, D4d, D5d, D6d, D7d, D8d,
  T, Td, Th,
  O, Oh,
  I, Ih,
  Cinfv, Dinfh
};

} // namespace Shapes
} // namespace Molassembler
} // namespace Scine

#endif
