#ifndef INCLUDE_MOLECULE_IO_H
#define INCLUDE_MOLECULE_IO_H

#include "detail/FileHandlers.h"

/*! @file
 *
 * Contains main IO definitions of the library. Currently only supports
 * MOLFile V2000 specification.
 */

/* TODO
 * - test V2000
 * - implement V3000
 * - Read coordinates are currently discarded. Either process them for
 *   information or use them to check against the CTAB
 */

namespace molassembler {

//! Input and output
namespace IO {

/*! Read a single molecule from a file.
 *
 * @throws If interpretation of coordinates and connectivity yields multiple
 *   molecules.
 * @note Interprets file type from extension. mol is a MOLFile, xyz an XYZ file
 *   and masm a CBOR serial representation of Molecule
 */
Molecule read(const std::string& filename);

/*! Read multiple molecules from a file.
 *
 * @note Interprets file type from extension. mol is a MOLFile, xyz an XYZ file
 *   and masm a CBOR serial representation of Molecule
 */
std::vector<Molecule> split(const std::string& filename);

/*! Writer function for MOL and XYZ formats
 *
 * @throws If the file extension does not match either mol or xyz.
 * @note Interprets which file type is to be written from filename extension.
 */
void write(
  const std::string& filename,
  const Molecule& molecule,
  const AngstromWrapper& angstromWrapper,
  const FileHandlers::IndexPermutation permutation = FileHandlers::IndexPermutation::Identity
);

//! Writer function from a PositionCollection in bohr
void write(
  const std::string& filename,
  const Molecule& molecule,
  const Delib::PositionCollection& positions,
  const FileHandlers::IndexPermutation permutation = FileHandlers::IndexPermutation::Identity
);

/*! Writer function to make files containing binary Molecule representation.
 *
 * @throws If the file extension does not match .masm
 */
void write(
  const std::string& filename,
  const Molecule& molecule
);

} // namespace IO

} // namespace molassembler

#endif
