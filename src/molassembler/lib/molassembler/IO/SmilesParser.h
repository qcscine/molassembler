/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Parse SMILES strings into molecules
 */
#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_PARSER_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_PARSER_H

#include "molassembler/Export.h"
#include <string>
#include <vector>

namespace Scine {
namespace molassembler {

class Molecule;

namespace io {
namespace experimental {

/**
 * @brief Parse a smiles string
 *
 * The smiles parser is implemented according to the OpenSMILES spec. It
 * supports the following features:
 * - Arbitrarily many molecules in a string
 * - Isotope markers (as long as they exist in Utils::ElementType)
 * - Valence filling of the organic subset
 * - Set shapes from VSEPR using (possibly) supplied charge
 * - Ring closures
 * - Stereo markers
 *   - Double bond
 *   - Tetrahedral (\@ / \@\@ / \@TH1 / \@TH2)
 *   - Square planar (\@SP1 - \@SP3)
 *   - Trigonal bipyramidal (\@TB1 - \@TB20)
 *   - Octahedral (\@OH1 - \@OH30)
 *
 * @warning Currently unsupported are: allene stereo markers and aromaticity.
 *
 * @param smiles the smiles string to parse containing possibly multiple
 *   molecules
 * @throws std::runtime_error If there are errors parsing the smiles string
 *
 * @return A list of molecules
 */
MASM_EXPORT std::vector<Molecule> parseSmiles(const std::string& smiles);

/**
 * @brief Parse a smiles string containing only a single molecule
 *
 * The smiles parser is implemented according to the OpenSMILES spec. It
 * supports the following features:
 * - Arbitrarily many molecules in a string
 * - Isotope markers (as long as they exist in Utils::ElementType)
 * - Valence filling of the organic subset
 * - Set shapes from VSEPR using (possibly) supplied charge
 * - Ring closures
 * - Stereo markers
 *   - Double bond
 *   - Tetrahedral (\@ / \@\@ / \@TH1 / \@TH2)
 *   - Square planar (\@SP1 - \@SP3)
 *   - Trigonal bipyramidal (\@TB1 - \@TB20)
 *   - Octahedral (\@OH1 - \@OH30)
 *
 * @warning Currently unsupported are: allene stereo markers and aromaticity.
 *
 * @param smiles The smiles string containing only a single molecule
 *
 * @throws std::logic_error If there are multiple molecules in the passed string
 * @throws std::runtime_error If there are errors parsing the smiles string
 *
 * @return A single molecule
 */
MASM_EXPORT Molecule parseSmilesSingleMolecule(const std::string& smiles);

} // namespace experimental
} // namespace io
} // namespace molassembler
} // namespace Scine

#endif
