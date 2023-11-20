/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
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

/**
 * @brief Generate smiles string for molecule
 *
 * @warning This is a lossy serialization format! The openSMILES standard does
 * not contain stereodescriptors for shapes other than the tetrahedron, square,
 * trigonal bipyramid and octahedron. Generated smiles containing stereocenters
 * with other shapes will not contain stereodescriptors for these centers.
 *
 * @note Missing normalization: Aromaticity detection in kekulized graph to aromatic
 * atom types.
 *
 * @returns A (partially) normalized openSMILES-standard compliant smiles
 *   string.
 */
MASM_EXPORT std::string emitSmiles(const Molecule& mol);

} // namespace Experimental
} // namespace IO
} // namespace Molassembler
} // namespace Scine

#endif
