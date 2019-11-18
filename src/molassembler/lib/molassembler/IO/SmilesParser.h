/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Parse SMILES strings into molecules
 */
#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_PARSER_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_PARSER_H

#include <string>
#include <vector>

namespace Scine {
namespace molassembler {

class Molecule;

namespace IO {

std::vector<Molecule> parseSmiles(const std::string& smiles);

} // namespace IO

} // namespace molassembler
} // namespace Scine

#endif
