/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/IO.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "molassembler/IO/BinaryHandler.h"
#include "molassembler/Interpret.h"
#include "molassembler/Molecule.h"
#include "molassembler/Options.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/Serialization.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"
#include "Utils/IO/ChemicalFileFormats/OpenBabelStreamHandler.h"

#include "temple/Random.h"

#include <fstream>
#include <iomanip>
#include <ctime>

namespace Scine {

namespace molassembler {

namespace IO {

const bool& LineNotation::enabled() {
  static bool enable = Utils::OpenBabelStreamHandler::checkForBinary();
  return enable;
}

Molecule LineNotation::fromFormat(const std::string& lineNotation, const std::string& format) {
  if(!enabled()) {
    throw std::logic_error(
      "obabel was not found! Line notation Molecule sources are unavailable."
    );
  }

  std::stringstream stream(lineNotation);
  Utils::OpenBabelStreamHandler handler;
  auto data = handler.read(stream, format);
  auto interpretation = interpret::molecules(
    data.first,
    data.second,
    interpret::BondDiscretizationOption::RoundToNearest
  );

  if(interpretation.molecules.size() > 1) {
    throw std::logic_error(
      "Interpreted multiple molecules into obabel output. Have you supplied multiple molecules in the line notation?"
    );
  }

  return interpretation.molecules.front();
}

Molecule LineNotation::fromCanonicalSMILES(const std::string& can) {
  return fromFormat(can, "can");
}

Molecule LineNotation::fromIsomericSMILES(const std::string& smi) {
  return fromFormat(smi, "smi");
}

Molecule LineNotation::fromInChI(const std::string& inchi) {
  return fromFormat(inchi, "inchi");
}

std::pair<Utils::AtomCollection, Utils::BondOrderCollection> exchangeFormat(
  const Molecule& molecule,
  AngstromWrapper angstromWrapper
) {
  return std::make_pair(
    Utils::AtomCollection(
      molecule.graph().elementCollection(),
      angstromWrapper.getBohr()
    ),
    molecule.graph().bondOrders()
  );
}

std::pair<Utils::AtomCollection, Utils::BondOrderCollection> exchangeFormat(
  const Molecule& molecule,
  const Utils::PositionCollection& positions
) {
  return std::make_pair(
    Utils::AtomCollection(
      molecule.graph().elementCollection(),
      positions
    ),
    molecule.graph().bondOrders()
  );
}

std::tuple<Utils::AtomCollection, Utils::BondOrderCollection, std::vector<AtomIndex>> shuffle(
  const Utils::AtomCollection& ac,
  const Utils::BondOrderCollection& bos
) {
  const unsigned N = ac.size();

  std::vector<AtomIndex> permutation;
  permutation.resize(N);
  std::iota(std::begin(permutation), std::end(permutation), 0);
  temple::random::shuffle(permutation, randomnessEngine());

  Utils::AtomCollection permutedAtoms(N);
  for(unsigned i = 0; i < N; ++i) {
    permutedAtoms.setPosition(permutation.at(i), ac.getPosition(i));
    permutedAtoms.setElement(permutation.at(i), ac.getElement(i));
  }

  Utils::BondOrderCollection permutedBOs(N);
  using SparseMatrixType = std::decay_t<
    decltype(std::declval<Utils::BondOrderCollection>().getMatrix())
  >;
  const SparseMatrixType& boMatrix = bos.getMatrix();
  // Iterate through the sparse representation
  for(int k = 0; k < boMatrix.outerSize(); ++k) {
    for(SparseMatrixType::InnerIterator it(boMatrix, k); it; ++it) {
      permutedBOs.setOrder(
        permutation.at(it.row()),
        permutation.at(it.col()),
        it.value()
      );
    }
  }

  return std::make_tuple(
    std::move(permutedAtoms),
    std::move(permutedBOs),
    permutation
  );
}

Molecule read(const std::string& filename) {
  boost::filesystem::path filepath {filename};
  if(!boost::filesystem::exists(filepath)) {
    throw std::logic_error("File selected to read does not exist.");
  }

  // Direct serializations of molecules have their own filetypes
  if(filepath.extension() == ".cbor") {
    return JsonSerialization(
      BinaryHandler::read(filename),
      JsonSerialization::BinaryFormat::CBOR
    );
  }

  if(filepath.extension() == ".bson") {
    return JsonSerialization(
      BinaryHandler::read(filename),
      JsonSerialization::BinaryFormat::BSON
    );
  }

  if(filepath.extension() == ".json") {
    std::ifstream input(filename);
    std::stringstream buffer;
    buffer << input.rdbuf();
    Molecule mol = JsonSerialization(buffer.str());
    input.close();
    return mol;
  }

  // This can throw in lots of cases
  auto readData = Utils::ChemicalFileHandler::read(filename);

  interpret::MoleculesResult interpretation;
  if(readData.second.empty()) {
    // Unfortunately, the file type does not include bond order information
    interpretation = interpret::molecules(readData.first, interpret::BondDiscretizationOption::RoundToNearest);
  } else {
    interpretation = interpret::molecules(readData.first, readData.second, interpret::BondDiscretizationOption::RoundToNearest);
  }

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

  // This can throw in lots of cases
  auto readData = Utils::ChemicalFileHandler::read(filename);

  interpret::MoleculesResult interpretation;
  if(readData.second.empty()) {
    // Unfortunately, the file type does not include bond order information
    interpretation = interpret::molecules(readData.first);
  } else {
    interpretation = interpret::molecules(readData.first, readData.second);
  }

  return interpretation.molecules;
}

void write(
  const std::string& filename,
  const Molecule& molecule,
  const AngstromWrapper& angstromWrapper
) {
  assert(molecule.graph().N() == static_cast<AtomIndex>(angstromWrapper.positions.rows()));
  auto data = exchangeFormat(molecule, angstromWrapper);
  Utils::ChemicalFileHandler::write(filename, data.first, data.second);
}

void write(
  const std::string& filename,
  const Molecule& molecule,
  const Utils::PositionCollection& positions
) {
  auto data = exchangeFormat(molecule, positions);

  // TODO Utils throws if bond information is supplied but the required format doesnt store it. Where to fix?
  Utils::ChemicalFileHandler::write(filename, data.first, data.second);
}

void write(
  const std::string& filename,
  const Molecule& molecule
) {
  boost::filesystem::path filepath {filename};

  if(filepath.extension() == ".cbor") {
    BinaryHandler::write(
      filename,
      JsonSerialization(molecule).toBinary(JsonSerialization::BinaryFormat::CBOR)
    );
    return;
  }

  if(filepath.extension() == ".bson") {
    BinaryHandler::write(
      filename,
      JsonSerialization(molecule).toBinary(JsonSerialization::BinaryFormat::BSON)
    );
    return;
  }

  if(filepath.extension() == ".json") {
    std::ofstream outfile(filename);
    outfile << JsonSerialization(molecule).operator std::string();
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
