/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_COMMON_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_COMMON_H

#include "Utils/Geometry/ElementInfo.h"
#include "Molassembler/Graph/PrivateGraph.h"

namespace Scine {
namespace Molassembler {
namespace IO {

bool isValenceFillElement(Utils::ElementType e);

int vertexValence(PrivateGraph::Vertex i, const PrivateGraph& g);

unsigned valenceFillElementImplicitHydrogenCount(
  int valence,
  Utils::ElementType e
);

} // namespace IO
} // namespace Molassembler
} // namespace Scine

#endif
