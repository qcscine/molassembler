#ifndef INCLUDE_MOLASSEMBLER_IO_H
#define INCLUDE_MOLASSEMBLER_IO_H

#include "molassembler/Types.h"

#include <string>
#include <vector>

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

// Forward declarations
namespace Delib {
class PositionCollection;
} // namespace Delib

namespace molassembler {

// Forward declarations
class Molecule;
class AngstromWrapper;

//! Input and output
namespace IO {

//! Defines permutations of atom indices of the molecule
enum class IndexPermutation {
  //! Atom indices are unaltered
  Identity,
  //! Atom indices are sorted by element type
  SortByElement,
  //! Randomize atom indices
  Random
};

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
  IndexPermutation permutation = IndexPermutation::Identity
);

//! Writer function from a PositionCollection in bohr
void write(
  const std::string& filename,
  const Molecule& molecule,
  const Delib::PositionCollection& positions,
  IndexPermutation permutation = IndexPermutation::Identity
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
