/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_COMMON_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_COMMON_H

#include "Utils/Geometry/ElementInfo.h"

namespace Scine {
namespace Molassembler {
namespace IO {

bool isValenceFillElement(Utils::ElementType e);

unsigned valenceFillElementImplicitHydrogenCount(
  int valence,
  Utils::ElementType e,
  bool aromatic
);

} // namespace IO
} // namespace Molassembler
} // namespace Scine

#endif
