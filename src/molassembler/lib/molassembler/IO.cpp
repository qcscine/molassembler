/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/IO.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "molassembler/OuterGraph.h"
#include "molassembler/IO/FileHandlers.h"
#include "molassembler/Molecule.h"
#include "molassembler/Serialization.h"

#include <fstream>
#include <iomanip>
#include <ctime>

namespace Scine {

namespace molassembler {

namespace IO {

Molecule read(const std::string& filename) {
  boost::filesystem::path filepath {filename};
  if(!boost::filesystem::exists(filepath)) {
    throw std::logic_error("File selected to read does not exist.");
  }

  if(filepath.extension() == ".masm") {
    return fromCBOR(
      FileHandlers::BinaryHandler::read(filename)
    );
  }

  if(filepath.extension() == ".json") {
    std::ifstream input(filename);
    std::stringstream buffer;
    buffer << input.rdbuf();
    auto mol = fromJSON(buffer.str());
    input.close();
    return mol;
  }

  std::unique_ptr<FileHandlers::FileHandler> handler;
  if(filepath.extension() == ".mol") {
    handler = std::make_unique<FileHandlers::MOLFileHandler>();
  } else if(filepath.extension() == ".xyz") {
    handler = std::make_unique<FileHandlers::XYZHandler>();
  } else {
    throw std::logic_error("Can only read files with .mol or .xyz extensions!");
  }

  auto interpretation = FileHandlers::interpret(
    handler->read(filename)
  );

  if(interpretation.molecules.size() > 1) {
    throw std::runtime_error(
      std::string("File is not a single molecule, but contains ")
        + std::to_string(interpretation.molecules.size())
        + " components."
    );
  }

  return interpretation.molecules.front();
}

std::vector<Molecule> split(const std::string& filename) {
  boost::filesystem::path filepath {filename};
  if(!boost::filesystem::exists(filepath)) {
    throw std::logic_error("File selected to read does not exist.");
  }

  std::unique_ptr<FileHandlers::FileHandler> handler;
  if(filepath.extension() == ".mol") {
    handler = std::make_unique<FileHandlers::MOLFileHandler>();
  } else if(filepath.extension() == ".xyz") {
    handler = std::make_unique<FileHandlers::XYZHandler>();
  } else {
    throw std::logic_error("Can only read files with .mol or .xyz extensions!");
  }

  auto interpretation = FileHandlers::interpret(
    handler->read(filename)
  );

  return interpretation.molecules;
}

void write(
  const std::string& filename,
  const Molecule& molecule,
  const AngstromWrapper& angstromWrapper,
  const IndexPermutation permutation
) {
  assert(molecule.graph().N() == static_cast<AtomIndex>(angstromWrapper.positions.size()));

  boost::filesystem::path filepath {filename};

  std::unique_ptr<FileHandlers::FileHandler> handler;
  if(filepath.extension() == ".mol") {
    handler = std::make_unique<FileHandlers::MOLFileHandler>();
  } else if(filepath.extension() == ".xyz") {
    handler = std::make_unique<FileHandlers::XYZHandler>();
  } else {
    throw std::logic_error("Can only read files with .mol or .xyz extensions!");
  }

  handler->write(filename, molecule, angstromWrapper, permutation);
}

void write(
  const std::string& filename,
  const Molecule& molecule,
  const Scine::Utils::PositionCollection& positions,
  const IndexPermutation permutation
) {
  AngstromWrapper wrapper {positions};
  return write(filename, molecule, wrapper, permutation);
}

void write(
  const std::string& filename,
  const Molecule& molecule
) {
  boost::filesystem::path filepath {filename};

  if(filepath.extension() == ".masm") {
    FileHandlers::BinaryHandler::write(
      filename,
      toCBOR(molecule)
    );
    return;
  }

  if(filepath.extension() == ".json") {
    std::ofstream outfile(filename);
    outfile << toJSON(molecule);
    outfile.close();
    return;
  }

  throw std::logic_error(
    "It makes no sense to write MOL or XYZ files without a PositionCollection"
  );
}

} // namespace IO

} // namespace molassembler

} // namespace Scine
