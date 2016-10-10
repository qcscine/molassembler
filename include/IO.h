#ifndef INCLUDE_MOLECULE_IO_H
#define INCLUDE_MOLECULE_IO_H

#include "Molecule.h"
#include "GraphAlgorithms.h"

#include <fstream>

/* TODO
 * - test V2000
 * - implement V3000
 */

namespace MoleculeManip {

namespace IO {

class FileReader {
public:
  virtual bool canReadFile(const std::string& filename) = 0;
  virtual Molecule readSingle(const std::string& filename) = 0;
};

class MOLFileReader : public FileReader {
private:
/* Typedefs */
  enum State {
    HeaderMoleculeName,
    HeaderProgram,
    HeaderComments,
    CountsLine,
    AtomBlock,
    BondBlock,
    End
  };

  enum MOLFileVersion {
    V2000
    // V3000
  };


/* Private members */
  Delib::ElementTypeCollection _elements;
  Delib::PositionCollection _positions;
  AdjacencyList _adjacencies;
  EdgeList _edges;

  const std::map<std::string, MOLFileVersion> _versionMap {
    {"V2000", MOLFileVersion::V2000}
    // {"V3000", MOLFileVersion::V3000}
  };

  const std::map<unsigned, BondType> _bondTypeMap {
    {1, BondType::Single},
    {2, BondType::Double},
    {3, BondType::Triple},
    {4, BondType::Aromatic}
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

      // Element name
      _elements.push_back(
        Delib::ElementInfo::instance()[
          _removeAllSpaces(
            line.substr(31, 3)
          )
        ].first
      );

      // Update adjacencies
      _adjacencies.addSlot();
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

      // WARNING: this can throw if queries are set!
      // update edges
      _edges.add(Edge(
        i,
        j,
        _bondTypeMap.at(bty)
      ));

      // add adjacencies
      _adjacencies.addAdjacency(i, j);
    }
  }

public:
/* Public member functions */
  virtual bool canReadFile(const std::string& filename) override {
    return (
      filename.substr(
        filename.length() - 4
      ) == ".mol"
    );
  }

  /*!
   * Throws in a myriad of cases!
   */
  virtual Molecule readSingle(const std::string& filename) override {
    assert(canReadFile(filename));

    std::ifstream file(filename);

    // Initialization 
    State state = State::HeaderMoleculeName;
    MOLFileVersion formatVersion = MOLFileVersion::V2000;
    // clear private state
    _positions.clear();
    _elements.clear();
    _adjacencies.clear();
    _edges.clear();
    
    unsigned atomBlockSize, bondBlockSize;
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
    unsigned nComponents = GraphAlgorithms::num_connected_components(
      _adjacencies
    );
    if(nComponents != 1) {
      throw std::runtime_error(
        std::string("File is not a single molecule, but contains ")
          + std::to_string(nComponents)
          + " components."
      );
    }

    return Molecule(
      _elements,
      _positions,
      _adjacencies,
      _edges
    );
  }
};

} // eo namespace IO

}

#endif
