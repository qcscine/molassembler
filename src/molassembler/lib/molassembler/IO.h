/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Input/output interface
 *
 * Contains main IO definitions of the library. Currently only supports
 * MOLFile V2000 specification.
 *
 * @todo
 * - test MOLFile V2000
 * - implement MOLFile V3000
 */

#ifndef INCLUDE_MOLASSEMBLER_IO_H
#define INCLUDE_MOLASSEMBLER_IO_H

#include "molassembler/Types.h"
#include "Utils/Typenames.h"

#include <string>
#include <vector>

namespace Scine {

// Forward declarations
namespace Utils {
class AtomCollection;
class BondOrderCollection;
} // namespace Utils

namespace molassembler {

// More forward declarations
class Molecule;
class AngstromWrapper;

//! Input and output
namespace IO {

/**
 * @brief Provides Molecule instances from line notations of molecules such
 *   as SMILES and InChI
 */
class LineNotation {
public:
  //! Checks whether the `obabel` binary is found in your path.
  static const bool& enabled();
  //! Construct a single molecule from a canonical SMILES string
  static Molecule fromCanonicalSMILES(const std::string& can);
  //! Construct a single molecule from an isomeric SMILES string
  static Molecule fromIsomericSMILES(const std::string& smi);
  /*!
   * @brief Construct a single molecule from an InChI string
   * @note The passed string has to include the `InChI=` prefix
   */
  static Molecule fromInChI(const std::string& inchi);

private:
  static Molecule fromFormat(const std::string& lineNotation, const std::string& format);
};

//! Extract exchange format information from a molecule and positional data
std::pair<Utils::AtomCollection, Utils::BondOrderCollection> exchangeFormat(
  const Molecule& molecule,
  AngstromWrapper angstromWrapper
);

//! @overload
std::pair<Utils::AtomCollection, Utils::BondOrderCollection> exchangeFormat(
  const Molecule& molecule,
  const Utils::PositionCollection& positions
);

//! Applies a random atom index permutation to exchange data
std::pair<Utils::AtomCollection, Utils::BondOrderCollection> shuffle(
  const Utils::AtomCollection& ac,
  const Utils::BondOrderCollection& bos
);

/*!
 * @brief Read a single molecule from a file.
 * @throws If interpretation of coordinates and connectivity yields multiple
 *   molecules.
 * @note Interprets file type from extension. mol is a MOLFile, xyz an XYZ file
 *   and masm a CBOR serial representation of Molecule
 */
Molecule read(const std::string& filename);

/*!
 * @brief Read multiple molecules from a file.
 *
 * @note Interprets file format from its extension.
 * @note masm and json serializations of Molecules cannot be split, they always
 *   contain only a single molecule. Use @p read() instead.
 */
std::vector<Molecule> split(const std::string& filename);

/*!
 * @brief Writer function for various chemical formats
 *
 * For exceptions this might throw, @see Utils::ChemicalFileHandler::write
 *
 * @note Interprets which file type is to be written from filename extension.
 */
void write(
  const std::string& filename,
  const Molecule& molecule,
  const AngstromWrapper& angstromWrapper
);

//! @overload
void write(
  const std::string& filename,
  const Molecule& molecule,
  const Utils::PositionCollection& positions
);

/*!
 * @brief Writer function for Molecule serializations
 * @throws If the file extension does not match .masm or .json
 */
void write(
  const std::string& filename,
  const Molecule& molecule
);

} // namespace IO

} // namespace molassembler

} // namespace Scine

#endif
