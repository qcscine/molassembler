/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/IO.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/process/child.hpp"
#include "boost/process/io.hpp"
#include "boost/process/search_path.hpp"

#include "Molassembler/IO/BinaryHandler.h"
#include "Molassembler/Interpret.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Options.h"
#include "Molassembler/Graph.h"
#include "Molassembler/Serialization.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"
#include "Utils/IO/ChemicalFileFormats/OpenBabelStreamHandler.h"

#include "Molassembler/Temple/Random.h"

#include <fstream>
#include <iomanip>
#include <ctime>

namespace Scine {
namespace Molassembler {
namespace IO {
namespace {

std::string pipeSvg(const Scine::Molassembler::Molecule& m) {
  std::stringstream os;

  // Construct pipe streams for redirection
  boost::process::opstream ips;
  boost::process::pstream ps;
  boost::process::pstream err;

  // Start the child process
  boost::process::child childProcess("dot -Tsvg", boost::process::std_in<ips, boost::process::std_out> ps,
                                     boost::process::std_err > err);

  // Feed our graphviz into the process
  ips << m.dumpGraphviz();
  ips.flush();
  ips.pipe().close();

  // Wait for the child process to exit
  childProcess.wait();

  std::stringstream stderrStream;
#if BOOST_VERSION >= 107000
  /* NOTE: This implementation of buffer transfers in boost process has a bug
   * that isn't fixed before Boost 1.70.
   */
  os << ps.rdbuf();
  stderrStream << err.rdbuf();
#else
  // Workaround: cast to a parent class implementing rdbuf() correctly.
  using BasicIOSReference = std::basic_ios<char, std::char_traits<char>>&;
  // Feed the results into our ostream
  os << static_cast<BasicIOSReference>(ps).rdbuf();
  stderrStream << static_cast<BasicIOSReference>(err).rdbuf();
#endif

  return os.str();
}

} // namespace

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
  auto interpretation = Interpret::molecules(
    data.first,
    data.second,
    Interpret::BondDiscretizationOption::RoundToNearest
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
  const AngstromPositions& angstromWrapper
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
  Temple::Random::shuffle(permutation, randomnessEngine());

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

  Interpret::MoleculesResult interpretation;
  if(readData.second.empty()) {
    // Unfortunately, the file type does not include bond order information
    interpretation = Interpret::molecules(readData.first, Interpret::BondDiscretizationOption::RoundToNearest);
  } else {
    interpretation = Interpret::molecules(readData.first, readData.second, Interpret::BondDiscretizationOption::RoundToNearest);
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

  Interpret::MoleculesResult interpretation;
  if(readData.second.empty()) {
    // Unfortunately, the file type does not include bond order information
    interpretation = Interpret::molecules(readData.first);
  } else {
    interpretation = Interpret::molecules(readData.first, readData.second);
  }

  return interpretation.molecules;
}

void write(
  const std::string& filename,
  const Molecule& molecule,
  const AngstromPositions& angstromWrapper
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

void write(const std::string& filename, const Molecule& molecule) {
  boost::filesystem::path filepath {filename};

  if(filepath.extension() == ".dot") {
    std::ofstream dotfile(filename);
    dotfile << molecule.dumpGraphviz();
    dotfile.close();
    return;
  }

  if(filepath.extension() == ".svg") {
    if(boost::process::search_path("dot").empty()) {
      throw std::runtime_error("Graphviz 'dot' binary not found in PATH");
    } else {
      std::ofstream svgfile(filename);
      svgfile << pipeSvg(molecule);
      svgfile.close();
      return;
    }
  }

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
    "File suffix cannot be mapped to a valid filetype for this writer"
  );
}

} // namespace IO
} // namespace Molassembler
} // namespace Scine
