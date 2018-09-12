#include "molassembler/IO/FileHandlers.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/range/iterator_range_core.hpp"

#include <fstream>
#include <iomanip>
#include <ctime>

#include "Delib/AtomCollection.h"
#include "Delib/AtomCollectionIO.h"
#include "Delib/Constants.h"
#include "Delib/ElementInfo.h"

#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Options.h"
#include "molassembler/Version.h"

namespace molassembler {

namespace IO {

namespace permutations {

Random::Random(AtomIndex N) {
  permutation.resize(N);
  std::iota(permutation.begin(), permutation.end(), 0);

  prng.shuffle(permutation);
}

SortByElement::SortByElement(const Molecule& mol) {
  permutation.resize(mol.graph().N());
  std::iota(permutation.begin(), permutation.end(), 0);
  std::sort(
    permutation.begin(),
    permutation.end(),
    [&mol](const AtomIndex i, const AtomIndex j) -> bool {
      return mol.graph().elementType(i) < mol.graph().elementType(j);
    }
  );
}

} // namespace permutations

namespace FileHandlers {

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

/* CTFile specification for V2000
 * 1 - Single
 * 2 - Double
 * 3 - Triple
 * 4 - Aromatic
 * 5 - Single / Double
 * 6 - Single / Aromatic
 * 7 - Double / Aromatic
 * 8 - Any
 */
const boost::bimap<unsigned, BondType> MOLFileHandler::_bondTypeMap
= makeBimap<unsigned, BondType>({
  {1, BondType::Single},
  {2, BondType::Double},
  {3, BondType::Triple},
  {8, BondType::Eta}
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
  const AngstromWrapper& angstromWrapper,
  const MOLFileVersion& version,
  const IndexPermutation& permutation
) const {
  std::function<AtomIndex(AtomIndex)> indexMap, inverse;

  const unsigned N = molecule.graph().N();

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
        << std::setw(8) << "MASM"+version::majorMinor() // PPPPPPPP (8, Prog name)
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
        << std::setw(3) << molecule.graph().B() // bbb
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
      auto symbolStringLambda = [](const Delib::ElementType elementType) {
        return Delib::ElementInfo::symbol(elementType);
      };

      for(unsigned i = 0; i < molecule.graph().N(); i++) {
        fout << std::setprecision(4) << std::fixed
          << std::setw(10) << angstromWrapper.positions[inverse(i)].x()
          << std::setw(10) << angstromWrapper.positions[inverse(i)].y()
          << std::setw(10) << angstromWrapper.positions[inverse(i)].z()
          << " " << std::setprecision(0)
          // aaa (atom symbol)
          << std::setw(3) << symbolStringLambda(molecule.graph().elementType(inverse(i)))
          << std::setw(2) << 0u // dd (isotope mass difference)
          << std::setw(3) << 0u // ccc (local charge) -> distant TODO
          << std::setw(3) << 0u // sss (atom stereo parity, ignored)
          << std::setw(3) << 0u // hhh (hydrogen count, for query, ignored)
          << std::setw(3) << 0u // bbb (stereo care box??, ignored)
          // vvv (valence)
          << std::setw(3) << molecule.graph().inner().degree(inverse(i))
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
      const InnerGraph& inner = molecule.graph().inner();

      for(
        const InnerGraph::Edge& edge :
        boost::make_iterator_range(
          inner.edges()
        )
      ) {
        auto i = indexMap(inner.source(edge));
        auto j = indexMap(inner.target(edge));
        auto bty = _bondTypeMap.right.at(inner.bondType(edge));

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

        data.elements.push_back(
          Delib::ElementInfo::elementTypeForSymbol(
            _removeAllSpaces(
              line.substr(31, 3)
            )
          )
        );

        data.angstromWrapper.positions.push_back(atomPosition);
      }

      // decrement remaning size
      --atomBlockSize;
      if(atomBlockSize == 0) {
        state = State::BondBlock;
      }
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
      if(bondBlockSize == 0) {
        state = State::End;
      }
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
  const AngstromWrapper& angstromWrapper,
  const IndexPermutation& permutation
) const {
  _write(
    filename,
    molecule,
    angstromWrapper,
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
  Delib::AtomCollection atoms = Delib::AtomCollectionIO::read(filename);

  data.elements = atoms.getElements();
  data.angstromWrapper = AngstromWrapper {
    atoms.getPositions(),
    LengthUnit::Bohr
  };

  return data;
}

void XYZHandler::write(
  const std::string& filename,
  const Molecule& molecule,
  const AngstromWrapper& angstromWrapper,
  const IndexPermutation& permutation
) const {
  std::function<AtomIndex(AtomIndex)> indexMap;

  const unsigned N = molecule.graph().N();

  if(permutation == IndexPermutation::Identity) {
    indexMap = permutations::Identity {};
  } else if(permutation == IndexPermutation::SortByElement) {
    indexMap = permutations::SortByElement {molecule};
  } else if(permutation == IndexPermutation::Random) {
    indexMap = permutations::Random {N};
  }

  Delib::AtomCollection ac;

  // Delib AtomCollection IO expects positions in bohr
  for(unsigned i = 0; i < N; ++i) {
    ac.push_back(
      Delib::Atom {
        molecule.graph().elementType(indexMap(i)),
        angstromWrapper.positions.at(indexMap(i)) * Delib::bohr_per_angstrom
      }
    );
  }

  Delib::AtomCollectionIO::write(filename, ac);
}

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

InterpretResult interpret(const FileHandler::RawData& data) {
  /* Some readers may not set bond orders at all. In those cases, bondOrders
   * has size zero.
   */
  if(data.bondOrders.getSystemSize() > 0) {
    return interpret(
      data.elements,
      data.angstromWrapper,
      data.bondOrders,
      BondDiscretizationOption::RoundToNearest
    );
  }

  return interpret(
    data.elements,
    data.angstromWrapper,
    BondDiscretizationOption::RoundToNearest
  );
}

} // namespace FileHandlers

} // namespace IO

} // namespace molassembler
