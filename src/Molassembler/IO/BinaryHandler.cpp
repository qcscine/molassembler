/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/IO/BinaryHandler.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

namespace Scine {
namespace Molassembler {
namespace IO {

bool BinaryHandler::canRead(const std::string& filename) {
  boost::filesystem::path filepath {filename};

  return (
    boost::filesystem::exists(filepath)
    && filepath.extension() == ".masm"
  );
}

void BinaryHandler::write(
  const std::string& filename,
  const BinaryType& binary
) {
  std::ofstream file(filename, std::ios::binary);

  std::size_t nElements = binary.size();

  file.write(reinterpret_cast<const char*>(&nElements), sizeof(nElements));
  file.write(reinterpret_cast<const char*>(&binary[0]), binary.size() * sizeof(std::uint8_t));

  file.close(); // Write EOF and close handle
}

BinaryHandler::BinaryType BinaryHandler::read(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);

  BinaryType data;

  std::size_t binarySize;
  file.read(reinterpret_cast<char*>(&binarySize), sizeof(binarySize));

  if(binarySize > 0) {
    data.resize(binarySize);
    file.read(reinterpret_cast<char*>(&data[0]), binarySize * sizeof(std::uint8_t));
  }

  file.close();

  return data;
}

} // namespace IO
} // namespace Molassembler
} // namespace Scine
