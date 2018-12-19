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

#include "Utils/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileHandler.h"
#include "Utils/IO/OpenBabelStreamHandler.h"

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
  InterpretResult interpretation = interpret(
    data.first,
    data.second,
    BondDiscretizationOption::RoundToNearest
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

std::pair<Utils::AtomCollection, Utils::BondOrderCollection> shuffle(
  const Utils::AtomCollection& ac,
  const Utils::BondOrderCollection& bos
) {
  const unsigned N = ac.size();

  std::vector<unsigned> permutation;
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
  const SparseMatrixType& BOMatrix = bos.getMatrix();
  // Iterate through the sparse representation
  for(int k = 0; k < BOMatrix.outerSize(); ++k) {
    for(SparseMatrixType::InnerIterator it(BOMatrix, k); it; ++it) {
      permutedBOs.setOrder(
        permutation.at(it.row()),
        permutation.at(it.col()),
        it.value()
      );
    }
  }

  return std::make_pair(
    std::move(permutedAtoms),
    std::move(permutedBOs)
  );
}

Molecule read(const std::string& filename) {
  boost::filesystem::path filepath {filename};
  if(!boost::filesystem::exists(filepath)) {
    throw std::logic_error("File selected to read does not exist.");
  }

  // Direct serializations of molecules have their own filetypes
  if(filepath.extension() == ".masm") {
    return fromCBOR(
      BinaryHandler::read(filename)
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

  // This can throw in lots of cases
  auto readData = Utils::ChemicalFileHandler::read(filename);

  InterpretResult interpretation;
  if(readData.second.empty()) {
    // Unfortunately, the file type does not include bond order information
    interpretation = interpret(readData.first, BondDiscretizationOption::RoundToNearest);
  } else {
    interpretation = interpret(readData.first, readData.second, BondDiscretizationOption::RoundToNearest);
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

  InterpretResult interpretation;
  if(readData.second.empty()) {
    // Unfortunately, the file type does not include bond order information
    interpretation = interpret(readData.first);
  } else {
    interpretation = interpret(readData.first, readData.second);
  }

  return interpretation.molecules;
}

void write(
  const std::string& filename,
  const Molecule& molecule,
  const AngstromWrapper& angstromWrapper
) {
  assert(molecule.graph().N() == static_cast<AtomIndex>(angstromWrapper.positions.size()));
  auto data = exchangeFormat(molecule, angstromWrapper);
  Utils::ChemicalFileHandler::write(filename, data.first, data.second);
}

void write(
  const std::string& filename,
  const Molecule& molecule,
  const Scine::Utils::PositionCollection& positions
) {
  auto data = exchangeFormat(molecule, positions);
  Utils::ChemicalFileHandler::write(filename, data.first, data.second);
}

void write(
  const std::string& filename,
  const Molecule& molecule
) {
  boost::filesystem::path filepath {filename};

  if(filepath.extension() == ".masm") {
    BinaryHandler::write(
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
