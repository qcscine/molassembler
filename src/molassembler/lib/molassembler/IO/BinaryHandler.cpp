/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/IO/BinaryHandler.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/range/iterator_range_core.hpp"

namespace Scine {

namespace molassembler {

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

  unsigned nElements = binary.size();
  write(file, nElements);

  for(const auto& element: binary) {
    write(file, element);
  }

  file.close(); // Write EOF and close handle
}

BinaryHandler::BinaryType BinaryHandler::read(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);

  BinaryType data;

  auto binarySize = read<unsigned>(file);

  if(binarySize > 0) {
    data.resize(binarySize);
    for(unsigned i = 0; i < binarySize; ++i) {
      data.at(i) = read<std::uint8_t>(file);
    }
  }

  file.close();

  return data;
}

} // namespace IO

} // namespace molassembler

} // namespace Scine
