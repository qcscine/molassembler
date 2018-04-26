#include "IO.h"
#include "Version.h"
#include "BondDistance.h"

#include "Delib/AtomCollectionIO.h"
#include "Delib/Constants.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include <fstream>
#include <iomanip>
#include <ctime>

namespace molassembler {

namespace IO {

/* MOLFileHandler implementation */
template <typename L, typename R>
boost::bimap<L, R> makeBimap(
  std::initializer_list<typename boost::bimap<L, R>::value_type> list
) {
  return boost::bimap<L, R>(list.begin(), list.end());
}

const boost::bimap<std::string, MOLFileHandler::MOLFileVersion> MOLFileHandler::_versionMap
= makeBimap<std::string, MOLFileHandler::MOLFileVersion>({
  {"V2000", MOLFileVersion::V2000}
  // {"V3000", MOLFileVersion::V3000}
});

const boost::bimap<unsigned, BondType> MOLFileHandler::_bondTypeMap
= makeBimap<unsigned, BondType>({
  {1, BondType::Single},
  {2, BondType::Double},
  {3, BondType::Triple},
  {4, BondType::Aromatic}
});

std::string MOLFileHandler::_removeAllSpaces(std::string&& a) {
  a.erase(
    std::remove(
      a.begin(),
      a.end(),
      ' '
    ),
    a.end()
  );

  return a;
}

std::string MOLFileHandler::_removeAllSpaces(const std::string& a) {
  std::string stringCopy = a;
  return _removeAllSpaces(std::move(stringCopy));
}

void MOLFileHandler::_write(
  const std::string& filename,
  const Molecule& molecule,
  const Delib::PositionCollection& positions,
  const MOLFileVersion& version,
  const IndexPermutation& permutation
) const {
  std::function<AtomIndexType(const AtomIndexType&)> indexMap, inverse;

  unsigned N = molecule.numAtoms();

  if(permutation == IndexPermutation::Identity) {
    indexMap = permutations::Identity {};
    inverse = permutations::Identity {};
  } else if(permutation == IndexPermutation::SortByElement) {
    auto forwardPermutation = permutations::SortByElement {molecule};
    inverse = permutations::Inverse {forwardPermutation, N};
    indexMap = std::move(forwardPermutation);
  } else if(permutation == IndexPermutation::Random) {
    auto forwardPermutation = permutations::Random {N};
    inverse = permutations::Inverse {forwardPermutation, N};
    indexMap = std::move(forwardPermutation);
  }

  std::ofstream fout(filename);
  fout << std::setprecision(0);

  State state = State::HeaderMoleculeName;

  while(state != State::End) {
    if(state == State::HeaderMoleculeName) {
      fout << "Unnamed Molecule" << std::endl;
      state = State::HeaderProgram;
    } else if(state == State::HeaderProgram) {
      auto now = time(nullptr);
      auto localNow = *localtime(&now);

      fout << std::setw(2) << "##" // First and last initial of user
        << std::setw(8) << "MASM"+Version::majorMinor() // PPPPPPPP (8, Prog name)
        << std::put_time(&localNow, "%m%d%y%H%M") // MMDDYYHHmm
        << "3D" // dd (dimensionality)
        // Missing:
        // SS (integer scaling factor)
        // ss (float scaling factor, 10 digits long: bbbb.aaaaa)
        // EE (energy, 12 digits long: sbbbbb.aaaaa)
        // RRRRRR (registry number)
        << std::endl;
      state = State::HeaderComments;
    } else if(state == State::HeaderComments) {
      fout << std::endl;
      state = State::CountsLine;
    } else if(state == State::CountsLine) {
      fout << std::setw(3) << N // aaa
        << std::setw(3) << molecule.numBonds() // bbb
        << std::setw(3) << 0u // lll (number of atom lists)
        << std::setw(3) << 0u // fff (obsolete)
        << std::setw(3) << 0u // ccc (chiral or not?) unhandled -> distant TODO
        << std::setw(3) << 0u // sss (num s-text entries, irrelevant here)
        << std::setw(3) << 0u // xxx (obsolete)
        << std::setw(3) << 0u // rrr (obsolete)
        << std::setw(3) << 0u // ppp (obsolete)
        << std::setw(3) << 0u // iii (obsolete)
        << std::setw(3) << 999u // mmm (num add. prop.s, unsupported, default 999)
        << std::setw(6) << _versionMap.right.at(version) // vvvvvv (Version string)
        << std::endl;
      state = State::AtomBlock;
    } else if(state == State::AtomBlock) {
      auto symbolStringLambda = [](const Delib::ElementType& elementType) {
        return Delib::ElementInfo::symbol(elementType);
      };

      for(unsigned i = 0; i < molecule.numAtoms(); i++) {
        fout << std::setprecision(4) << std::fixed
          << std::setw(10) << positions[inverse(i)].x()
          << std::setw(10) << positions[inverse(i)].y()
          << std::setw(10) << positions[inverse(i)].z()
          << " " << std::setprecision(0)
          // aaa (atom symbol)
          << std::setw(3) << symbolStringLambda(molecule.getElementType(inverse(i)))
          << std::setw(2) << 0u // dd (isotope mass difference)
          << std::setw(3) << 0u // ccc (local charge) -> distant TODO
          << std::setw(3) << 0u // sss (atom stereo parity, ignored)
          << std::setw(3) << 0u // hhh (hydrogen count, for query, ignored)
          << std::setw(3) << 0u // bbb (stereo care box??, ignored)
          // vvv (valence)
          << std::setw(3) << molecule.getNumAdjacencies(inverse(i))
          << std::setw(3) << 0u // HHH (H0 designator, ISIS/Desktop, ignored)
          << std::setw(3) << 0u // rrr (unused)
          << std::setw(3) << 0u // iii (unused)
          // mmm (atom-atom mapping number, for reactions, ignored)
          << std::setw(3) << 0u
          // nnn (inversion/retention flag, for reactions, ignored)
          << std::setw(3) << 0u
          // eee (exact change flag, for reactions, ignored)
          << std::setw(3) << 0u
          << std::endl;
      }
      state = State::BondBlock;
    } else if(state == State::BondBlock) {
      for(const auto& edge: molecule.iterateEdges()) {
        auto edgeVertices = molecule.vertices(edge);
        auto i = indexMap(edgeVertices.front());
        auto j = indexMap(edgeVertices.back());
        auto bty = _bondTypeMap.right.at(molecule.getBondType(edge));

        // MOLFile indices are 1-based, have to add 1 to internal indices!
        fout << std::setw(3) << (1 + std::min(i, j)) // 111 (index of 1st atom)
          << std::setw(3) << (1 + std::max(i, j))  // 222 (index of 2nd atom)
          << std::setw(3) << bty // ttt (bond type)
          << std::setw(3) << 0u // sss (bond stereo, ignored for now)
          << std::setw(3) << 0u // xxx (unused)
          << std::setw(3) << 0u // rrr (bond topology) -> distant TODO
          << std::setw(3) << 0u // ccc (reacting center status, ignored)
          << std::endl;
      }

      // final line
      fout << "M END";
      state = State::End;
    }
  }

  fout.close();
}

/*!
 * Throws in a myriad of cases!
 */
FileHandler::RawData MOLFileHandler::read(const std::string& filename) const {
  assert(canRead(filename));

  std::ifstream file(filename);

  // Initialization
  State state = State::HeaderMoleculeName;
  MOLFileVersion formatVersion = MOLFileVersion::V2000;

  RawData data;

  unsigned atomBlockSize = 0, bondBlockSize = 0;
  std::string line;

  while(state != State::End) {
    if(state == State::HeaderMoleculeName) {
      std::getline(file, line);
      state = State::HeaderProgram;
    } else if(state == State::HeaderProgram) {
      std::getline(file, line);
      state = State::HeaderComments;
    } else if(state == State::HeaderComments) {
      std::getline(file, line);
      state = State::CountsLine;
    } else if(state == State::CountsLine) {
      std::getline(file, line);

      atomBlockSize = std::atoi(line.substr(0, 3).c_str());
      bondBlockSize = std::atoi(line.substr(3, 3).c_str());
      formatVersion = _versionMap.left.at(
        _removeAllSpaces(
          line.substr(33)
        )
      );

      // Reserve space in data members
      data.bondOrders.resize(atomBlockSize);
      data.bondOrders.setZero();

      state = State::AtomBlock;
    } else if(state == State::AtomBlock) {
      std::getline(file, line);

      if(formatVersion == MOLFileVersion::V2000) {
        // Extract position
        Delib::Position atomPosition;
        atomPosition.x() = std::stod(line.substr(0, 10));
        atomPosition.y() = std::stod(line.substr(10, 10));
        atomPosition.z() = std::stod(line.substr(20, 10));

        data.atoms.push_back(
          Delib::Atom {
            Delib::ElementInfo::elementTypeForSymbol(
              _removeAllSpaces(
                line.substr(31, 3)
              )
            ),
            atomPosition
          }
        );
      }

      // decrement remaning size
      --atomBlockSize;
      if(atomBlockSize == 0) state = State::BondBlock;
    } else if(state == State::BondBlock) {
      std::getline(file, line);

      if(formatVersion == MOLFileVersion::V2000) {
        // Extract bond information
        unsigned i, j, bty;
        // MOLFile indices are 1-based, thus subtract one
        i = std::atoi(line.substr(0, 3).c_str()) - 1;
        j = std::atoi(line.substr(3, 3).c_str()) - 1;
        bty = std::atoi(line.substr(6, 3).c_str());

        /* The bond order must be translated from the integer given
         * to the internal enum type, and then to a double representation
         * for Molecule::interpret to properly handle.
         */
        data.bondOrders.setOrder(
          i,
          j,
          Bond::bondOrderMap.at(
            static_cast<unsigned>(
              _bondTypeMap.left.at(bty)
            )
          )
        );
      }

      // decrement
      --bondBlockSize;
      if(bondBlockSize == 0) state = State::End;
    }
  }

  return data;
}

bool MOLFileHandler::canRead(const std::string& filename) const {
  boost::filesystem::path filepath {filename};

  return (
    boost::filesystem::exists(filepath)
    && filepath.extension() == ".mol"
  );
}
void MOLFileHandler::write(
  const std::string& filename,
  const Molecule& molecule,
  const Delib::PositionCollection& positions,
  const IndexPermutation& permutation
) const {
  _write(
    filename,
    molecule,
    positions,
    MOLFileVersion::V2000,
    permutation
  );
}

/* XYZHandler implementation */
bool XYZHandler::canRead(const std::string& filename) const {
  boost::filesystem::path filepath {filename};

  return (
    boost::filesystem::exists(filepath)
    && filepath.extension() == ".xyz"
  );
}

FileHandler::RawData XYZHandler::read(const std::string& filename) const {
  assert(canRead(filename));

  FileHandler::RawData data;
  data.atoms = Delib::AtomCollectionIO::read(filename);
  // Delib IO makes everything bohr units, I need Angstroms
  data.atoms.setPositions(
    data.atoms.getPositions() * Delib::angstrom_per_bohr
  );

  return data;
}

void XYZHandler::write(
  const std::string& filename,
  const Molecule& molecule,
  const Delib::PositionCollection& positions,
  const IndexPermutation& permutation
) const {
  std::function<AtomIndexType(const AtomIndexType&)> indexMap;

  unsigned N = molecule.numAtoms();

  if(permutation == IndexPermutation::Identity) {
    indexMap = permutations::Identity {};
  } else if(permutation == IndexPermutation::SortByElement) {
    indexMap = permutations::SortByElement {molecule};
  } else if(permutation == IndexPermutation::Random) {
    indexMap = permutations::Random {N};
  }

  Delib::AtomCollection ac;

  for(unsigned i = 0; i < N; ++i) {
    ac.push_back(
      Delib::Atom {
        molecule.getElementType(indexMap(i)),
        Delib::Position {positions.at(indexMap(i))}
      }
    );
  }

  Delib::AtomCollectionIO::write(filename, ac);
}

namespace detail {

Molecule::InterpretResult interpret(const FileHandler::RawData& data) {
  /* Some readers may not set bond orders at all. In those cases, bondOrders
   * has size zero.
   */
  if(data.bondOrders.getSystemSize() > 0) {
    return Molecule::interpret(
      data.atoms,
      data.bondOrders,
      Molecule::BondDiscretizationOption::UFF
    );
  }

  return Molecule::interpret(
    data.atoms,
    Molecule::BondDiscretizationOption::UFF
  );
}

} // namespace detail

Molecule read(const std::string& filename) {
  boost::filesystem::path filepath {filename};
  if(!boost::filesystem::exists(filepath)) {
    throw std::logic_error("File selected to read does not exist.");
  }

  std::unique_ptr<FileHandler> handler;
  if(filepath.extension() == ".mol") {
    handler = std::make_unique<MOLFileHandler>();
  } else if(filepath.extension() == ".xyz") {
    handler = std::make_unique<XYZHandler>();
  } else {
    throw std::logic_error("Can only read files with .mol or .xyz extensions!");
  }

  auto interpretation = detail::interpret(
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

  std::unique_ptr<FileHandler> handler;
  if(filepath.extension() == ".mol") {
    handler = std::make_unique<MOLFileHandler>();
  } else if(filepath.extension() == ".xyz") {
    handler = std::make_unique<XYZHandler>();
  } else {
    throw std::logic_error("Can only read files with .mol or .xyz extensions!");
  }

  auto interpretation = detail::interpret(
    handler->read(filename)
  );

  return interpretation.molecules;
}

void write(
  const std::string& filename,
  const Molecule& molecule,
  const Delib::PositionCollection& positions,
  const FileHandler::IndexPermutation& permutation
) {
  boost::filesystem::path filepath {filename};

  std::unique_ptr<FileHandler> handler;
  if(filepath.extension() == ".mol") {
    handler = std::make_unique<MOLFileHandler>();
  } else if(filepath.extension() == ".xyz") {
    handler = std::make_unique<XYZHandler>();
  } else {
    throw std::logic_error("Can only read files with .mol or .xyz extensions!");
  }

  handler->write(filename, molecule, positions, permutation);
}

} // namespace IO

} // namespace molassembler
