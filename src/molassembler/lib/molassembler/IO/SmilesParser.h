#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_PARSER_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_PARSER_H

#include <string>

namespace Scine {
namespace molassembler {

class Molecule;

namespace IO {

Molecule parseSMILES(const std::string& smiles);

} // namespace IO

} // namespace molassembler
} // namespace Scine

#endif
