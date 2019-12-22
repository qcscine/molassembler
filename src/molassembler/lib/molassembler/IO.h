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
namespace Utils {

// Forward declarations
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
class MASM_EXPORT LineNotation {
public:
  /*! @brief Checks whether the `obabel` binary is found in your path.
   *
   * @complexity{@math{\Theta(1)}}
   */
  static const bool& enabled();
  /*! @brief Construct a single molecule from a canonical SMILES string
   *
   * @complexity{@math{\Theta(N)} presumably, depends on OpenBabel}
   */
  static Molecule fromCanonicalSMILES(const std::string& can);
  /*! @brief Construct a single molecule from an isomeric SMILES string
   *
   * @complexity{@math{\Theta(N)} presumably, depends on OpenBabel}
   */
  static Molecule fromIsomericSMILES(const std::string& smi);
  /*! @brief Construct a single molecule from an InChI string
   *
   * @complexity{@math{\Theta(N)} presumably, depends on OpenBabel}
   * @note The passed string has to include the `InChI=` prefix
   */
  static Molecule fromInChI(const std::string& inchi);

private:
  static Molecule fromFormat(const std::string& lineNotation, const std::string& format);
};

/*! @brief Extract exchange format information from a molecule and positional data
 *
 * @complexity{@math{\Theta(N + B)}}
 */
MASM_EXPORT std::pair<Utils::AtomCollection, Utils::BondOrderCollection> exchangeFormat(
  const Molecule& molecule,
  AngstromWrapper angstromWrapper
);

//! @overload
MASM_EXPORT std::pair<Utils::AtomCollection, Utils::BondOrderCollection> exchangeFormat(
  const Molecule& molecule,
  const Utils::PositionCollection& positions
);

/*! @brief Applies a random atom index permutation to exchange data
 *
 * @complexity{@math{\Theta(N + B)}}
 */
MASM_EXPORT std::tuple<Utils::AtomCollection, Utils::BondOrderCollection, std::vector<AtomIndex>> shuffle(
  const Utils::AtomCollection& ac,
  const Utils::BondOrderCollection& bos
);

/*! @brief Read a single molecule from a file.
 *
 * @complexity{@math{\Theta(N)} typically}
 * @throws If interpretation of coordinates and connectivity yields multiple
 *   molecules.
 * @note Interprets file type from extension. mol is a MOLFile, xyz an XYZ file
 *   and cbor/bson/json are serializations of Molecule
 */
MASM_EXPORT Molecule read(const std::string& filename);

/*! @brief Read multiple molecules from a file.
 *
 * @complexity{@math{\Theta(N)} typically}
 * @note Interprets file format from its extension. See read()
 * @note Serializations of Molecules cannot be split, they always
 *   contain only a single molecule. Use @p read() instead.
 */
MASM_EXPORT std::vector<Molecule> split(const std::string& filename);

/*! @brief Writer function for various chemical formats
 *
 * For exceptions this might throw, @see Utils::ChemicalFileHandler::write
 *
 * @complexity{@math{\Theta(N)}}
 * @note Interprets which file type is to be written from filename extension.
 *
 * @note Unless the file format is .json or .masm, canonicalization state of
 * Molecule instances is lost.
 */
MASM_EXPORT void write(
  const std::string& filename,
  const Molecule& molecule,
  const AngstromWrapper& angstromWrapper
);

//! @overload
MASM_EXPORT void write(
  const std::string& filename,
  const Molecule& molecule,
  const Utils::PositionCollection& positions
);

/*! @brief Writer function for Molecule serializations
 *
 * @complexity{@math{\Theta(V + E + A + B)}}
 * @note Canonicalization state is retained using the molecule serializations.
 * @throws If the file extension does not match .masm or .json
 */
MASM_EXPORT void write(
  const std::string& filename,
  const Molecule& molecule
);

} // namespace IO
} // namespace molassembler
} // namespace Scine

#endif
