#ifndef INCLUDE_MOLECULE_IO_H
#define INCLUDE_MOLECULE_IO_H

#include "Molecule.h"
#include "Version.h"

#include "temple/Random.h"

#include <fstream>
#include <iomanip>
#include <ctime>

/*! @file
 *
 * Contains main IO definitions of the library. Currently only supports
 * MOLFile V2000 specification.
 */

/* TODO
 * - test V2000
 * - implement V3000
 * - Read coordinates are currently discarded. Either process them for
 *   information or use them to check against the CTAB
 */

namespace molassembler {

//! Input and output
namespace IO {

namespace permutations {

struct Identity {
  template<typename T>
  T operator() (const T& t) const {
    return t;
  }
};

struct Random {
  std::vector<AtomIndexType> permutation;

  Random(AtomIndexType N) {
    permutation.resize(N);
    std::iota(permutation.begin(), permutation.end(), 0);
    std::shuffle(permutation.begin(), permutation.end(), temple::random.randomEngine);
  }

  AtomIndexType operator() (const AtomIndexType& i) const {
    return permutation.at(i);
  }
};

struct SortByElement {
  std::vector<AtomIndexType> permutation;

  SortByElement(const Molecule& mol) {
    permutation.resize(mol.numAtoms());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(
      permutation.begin(),
      permutation.end(),
      [&mol](const AtomIndexType& i, const AtomIndexType& j) -> bool {
        return mol.getElementType(i) < mol.getElementType(j);
      }
    );
  }

  AtomIndexType operator() (const AtomIndexType& i) const {
    return permutation.at(i);
  }
};

struct Inverse {
  std::vector<AtomIndexType> permutation;

  template<typename Permutation>
  Inverse(const Permutation& ante, const AtomIndexType& size) {
    permutation.resize(size);
    for(AtomIndexType i = 0; i < size; ++i) {
      permutation.at(ante(i)) = i;
    }
  }

  AtomIndexType operator() (const AtomIndexType& i) const {
    return permutation.at(i);
  }
};

} // namespace permutations

//! Abstract base class / interface for file reading classes
struct FileReader {
  virtual bool canReadFile(const std::string& filename) = 0;
  virtual Molecule readSingle(const std::string& filename) = 0;
};

//! Abstract base class / interface for file writing classes
struct FileWriter {
  enum class IndexPermutation {
    Identity,
    SortByElement,
    Random
  };

  virtual bool canWriteFile(const std::string& filename) = 0;
  virtual void writeSingle(
    const std::string& filename,
    const Molecule& molecule,
    const Delib::PositionCollection& positions,
    const IndexPermutation& permutation = IndexPermutation::Identity
  ) = 0;
};

/*!
 * Handles IO of MOLFiles
 */
class MOLFileHandler : public FileReader, public FileWriter {
private:
/* Typedefs */
  enum class State {
    HeaderMoleculeName,
    HeaderProgram,
    HeaderComments,
    CountsLine,
    AtomBlock,
    BondBlock,
    End
  };

  enum class MOLFileVersion {
    V2000
    // V3000
  };


/* Private members */
  Delib::PositionCollection _positions;
  GraphType _graph;

  const std::map<std::string, MOLFileVersion> _versionMap {
    {"V2000", MOLFileVersion::V2000}
    // {"V3000", MOLFileVersion::V3000}
  };

  const std::map<MOLFileVersion, std::string> _versionStringMap {
    {MOLFileVersion::V2000, "V2000"}
  };

  const std::map<unsigned, BondType> _bondTypeMap {
    {1, BondType::Single},
    {2, BondType::Double},
    {3, BondType::Triple},
    {4, BondType::Aromatic}
  };

  const std::map<BondType, unsigned> _bondTypeUnsignedMap {
    {BondType::Single, 1},
    {BondType::Double, 2},
    {BondType::Triple, 3},
    {BondType::Aromatic, 4}
  };

/* Private member functions */
  std::string _removeAllSpaces(const std::string& string) {
    std::string stringCopy = string;
    stringCopy.erase(
      std::remove(
        stringCopy.begin(),
        stringCopy.end(),
        ' '
      ),
      stringCopy.end()
    );
    return stringCopy;
  }

  void _parseAtomBlockLine(
    const MOLFileVersion& version,
    const std::string& line
  ) {
    if(version == MOLFileVersion::V2000) {
      // Extract position
      Delib::Position atomPosition;
      atomPosition.x() = std::stod(line.substr(0, 10));
      atomPosition.y() = std::stod(line.substr(10, 10));
      atomPosition.z() = std::stod(line.substr(20, 10));
      _positions.push_back(atomPosition);

      // Update adjacencies
      auto vertex = boost::add_vertex(_graph);
      _graph[vertex].elementType = Delib::ElementInfo::elementTypeForSymbol(
        _removeAllSpaces(
          line.substr(31, 3)
        )
      );
    }
  }

  void _parseBondBlockLine(
    const MOLFileVersion& version,
    const std::string& line
  ) {
    if(version == MOLFileVersion::V2000) {
      // Extract bond information
      unsigned i, j, bty;
      // MOLFile indices are 1-based, thus subtract one
      i = std::atoi(line.substr(0, 3).c_str()) - 1;
      j = std::atoi(line.substr(3, 3).c_str()) - 1;
      bty = std::atoi(line.substr(6, 3).c_str());

      // Add bond to graph
      auto edgeAddPair = boost::add_edge(i, j, _graph);
      _graph[edgeAddPair.first].bondType = _bondTypeMap.at(bty);
    }
  }

  /* NOTE: Important cases, perhaps to be altered
   * TODO If BondTypes of 4-6 are to be stored in a MOLFile, the program will
   * crash, no elegant failure exists.
   */
  void _writeSingle(
    const std::string& filename,
    const Molecule& molecule,
    const Delib::PositionCollection& positions,
    const MOLFileVersion& version,
    const IndexPermutation& permutation
  ) {
    assert(canWriteFile(filename));

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
          << std::setw(6) << _versionStringMap.at(version) // vvvvvv (Version string)
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
        for(const auto& edge: molecule.getEdges()) {
          auto i = indexMap(edge.first.first);
          auto j = indexMap(edge.first.second);
          auto bty = _bondTypeUnsignedMap.at(edge.second);

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

public:
/* Public member functions */
  /* Reading functions */
  bool canReadFile(const std::string& filename) final {
    return (
      filename.substr(
        filename.length() - 4
      ) == ".mol"
    );
  }

  /*!
   * Throws in a myriad of cases!
   */
  Molecule readSingle(const std::string& filename) final {
    assert(canReadFile(filename));

    std::ifstream file(filename);

    // Initialization
    State state = State::HeaderMoleculeName;
    MOLFileVersion formatVersion = MOLFileVersion::V2000;
    // clear private state
    _positions.clear();
    _graph.clear();

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
        formatVersion = _versionMap.at(
          _removeAllSpaces(
            line.substr(33)
          )
        );

        state = State::AtomBlock;
      } else if(state == State::AtomBlock) {
        std::getline(file, line);
       _parseAtomBlockLine(
          formatVersion,
          line
        );

        // decrement remaning size
        atomBlockSize--;
        if(atomBlockSize == 0) state = State::BondBlock;
      } else if(state == State::BondBlock) {
        std::getline(file, line);

        _parseBondBlockLine(
          formatVersion,
          line
        );

        // decrement
        bondBlockSize--;
        if(bondBlockSize == 0) state = State::End;
      }
    }

    // Ensure that the Molecule is connected, no fragments are contained
    unsigned nComponents = GraphAlgorithms::numConnectedComponents(_graph);

    if(nComponents != 1) {
      throw std::runtime_error(
        std::string("File is not a single molecule, but contains ")
          + std::to_string(nComponents)
          + " components."
      );
    }

    // Determine if Positions are dummies
    unsigned nZeroLengthLowerBound = 0;

    // One atom may be positioned on 0, 0, 0 (or close), but not two
    for(const auto& position : _positions) {
      if(position.asEigenVector().norm() <= 1e-14) {
        nZeroLengthLowerBound += 1;
        if(nZeroLengthLowerBound > 1) {
          break;
        }
      }
    }

    if(nZeroLengthLowerBound > 1) {
      /* Assume all are dummies, no 3D Information present -> no Stereocenter
       *  assignments
       */
      return Molecule(_graph);
    } else { // We have 3D information, try to infer stereocenters (if present)
      return Molecule(_graph, _positions);
    }
  }

  bool canWriteFile(const std::string& filename) final {
    return (
      filename.substr(
        filename.length() - 4
      ) == ".mol"
    );
  }

  void writeSingle(
    const std::string& filename,
    const Molecule& molecule,
    const Delib::PositionCollection& positions,
    const IndexPermutation& permutation = IndexPermutation::Identity
  ) final {
    _writeSingle(
      filename,
      molecule,
      positions,
      MOLFileVersion::V2000,
      permutation
    );
  }

  const Delib::PositionCollection& getPositionCollection() const {
    return _positions;
  }

  const GraphType& getGraph() const {
    return _graph;
  }
};

} // namespace IO

} // namespace molassembler

#endif
