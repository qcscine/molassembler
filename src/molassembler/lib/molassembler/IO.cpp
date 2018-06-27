#include "IO.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "Serialization.h"

#include <fstream>
#include <iomanip>
#include <ctime>

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
  const FileHandlers::IndexPermutation permutation
) {
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
  const Delib::PositionCollection& positions,
  const FileHandlers::IndexPermutation permutation
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
  } else {
    throw std::logic_error(
      "It makes no sense to write MOL or XYZ files without a PositionCollection"
    );
  }
}

} // namespace IO

} // namespace molassembler
