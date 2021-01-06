/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Emit SMILES from molecules
 */
#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_EMITTER_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_EMITTER_H

#include "Molassembler/Export.h"
#include <string>

namespace Scine {
namespace Molassembler {

class Molecule;

namespace IO {
namespace Experimental {

MASM_EXPORT std::string emitSmiles(const Molecule& mol);

} // namespace Experimental
} // namespace IO
} // namespace Molassembler
} // namespace Scine

#endif
