/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Binary Input/output
 */

#ifndef INCLUDE_MOLASSEMBLER_IO_BINARY_H
#define INCLUDE_MOLASSEMBLER_IO_BINARY_H

#include <fstream>
#include <string>
#include <vector>
#include <bitset>

namespace Scine {
namespace Molassembler {

// Forward-declarations
class Molecule;

namespace IO {

//! Binary file IO
struct BinaryHandler {
  using BinaryType = std::vector<std::uint8_t>;

  //! Checks whether the file extension is readable
  static bool canRead(const std::string& filename);

  /*! @brief Writes binary to a file
   *
   * @complexity{@math{\Theta(N)}}
   */
  static void write(
    const std::string& filename,
    const BinaryType& binary
  );

  /*! @brief Reads binary from a file
   *
   * @complexity{@math{\Theta(N)}}
   */
  static BinaryType read(const std::string& filename);
};

} // namespace IO
} // namespace Molassembler
} // namespace Scine

#endif
