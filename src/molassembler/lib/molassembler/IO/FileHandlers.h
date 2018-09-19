// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_FILE_HANDLERS_H
#define INCLUDE_MOLASSEMBLER_FILE_HANDLERS_H

#include "boost/bimap.hpp"

#include "Delib/ElementTypeCollection.h"
#include "Delib/BondOrderCollection.h"

#include "molassembler/AngstromWrapper.h"
#include "molassembler/IO.h"
#include "molassembler/Interpret.h"

/*!@file
 *
 * @brief Input/output of various file formats
 *
 * Centralizes input and output of several file formats
 * - MOL (V2000)
 * - XYZ
 * - MASM (Binary)
 */

namespace molassembler {

// Forward-declarations
class Molecule;

namespace IO {

namespace permutations {

//! Maps any type T provided to itself
struct Identity {
  template<typename T>
  T operator() (const T& t) const {
    return t;
  }
};

//! Randomizes atom indices.
struct Random {
  std::vector<AtomIndex> permutation;

  explicit Random(AtomIndex N);

  inline AtomIndex operator() (const AtomIndex i) const {
    return permutation.at(i);
  };
};

//! Sorts an element's indices by atom elements' Z
struct SortByElement {
  std::vector<AtomIndex> permutation;

  explicit SortByElement(const Molecule& mol);

  inline AtomIndex operator() (const AtomIndex i) const {
    return permutation.at(i);
  }
};

//! Provides the inverse permutation to any permutation
struct Inverse {
  std::vector<AtomIndex> permutation;

  template<typename Permutation>
  inline explicit Inverse(const Permutation& ante, const AtomIndex size) {
    permutation.resize(size);
    for(AtomIndex i = 0; i < size; ++i) {
      permutation.at(ante(i)) = i;
    }
  }

  inline AtomIndex operator() (const AtomIndex i) const {
    return permutation.at(i);
  }
};

} // namespace permutations

namespace FileHandlers {

//! Abstract base class to all file handlers
struct FileHandler {
  struct RawData {
    Delib::ElementTypeCollection elements;
    AngstromWrapper angstromWrapper;
    Delib::BondOrderCollection bondOrders;
  };

  // Virtualize destructor
  virtual ~FileHandler() = default;

  // Demand strong const guarantees for quality of implementation
  virtual bool canRead(const std::string& filename) const = 0;
  virtual RawData read(const std::string& filename) const = 0;
  virtual void write(
    const std::string& filename,
    const Molecule& molecule,
    const AngstromWrapper& angstromWrapper,
    const IndexPermutation& permutation = IndexPermutation::Identity
  ) const = 0;
};

//! MOL file IO
class MOLFileHandler final : public FileHandler {
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

  // Version string from file <-> internal enum capturing current version
  static const boost::bimap<std::string, MOLFileVersion> _versionMap;

  // MOLFile bond type value <-> internal bond type enum
  static const boost::bimap<unsigned, BondType> _bondTypeMap;

/* Private static functions */
  static std::string _removeAllSpaces(std::string&& a);
  static std::string _removeAllSpaces(const std::string& a);

/* Private member functions */

  /* NOTE: Important cases, perhaps to be altered
   * TODO If BondTypes of 4-6 are to be stored in a MOLFile, the program will
   * crash, no elegant failure exists.
   */
  void _write(
    const std::string& filename,
    const Molecule& molecule,
    const AngstromWrapper& angstromWrapper,
    const MOLFileVersion& version,
    const IndexPermutation& permutation
  ) const;

public:
/* Public types */
/* Public member functions */
  /* Reading functions */
  bool canRead(const std::string& filename) const final;

  RawData read(const std::string& filename) const final;

  //! Write a single Molecule to a file
  void write(
    const std::string& filename,
    const Molecule& molecule,
    const AngstromWrapper& angstromWrapper,
    const IndexPermutation& permutation = IndexPermutation::Identity
  ) const final;
};

//! XYZ file IO
struct XYZHandler final : public FileHandler {
  bool canRead(const std::string& filename) const final;

  RawData read(const std::string& filename) const final;

  //! Write a single Molecule to a file
  void write(
    const std::string& filename,
    const Molecule& molecule,
    const AngstromWrapper& angstromWrapper,
    const IndexPermutation& permutation = IndexPermutation::Identity
  ) const final;
};

//! Binary file IO
struct BinaryHandler {
  using BinaryType = std::vector<std::uint8_t>;

  template<typename T>
  static std::enable_if_t<
    std::is_unsigned<T>::value,
    void
  > write(
    std::ofstream& file,
    const T value
  ) {
    std::bitset<8 * sizeof(T)> bits {value};
    file << bits;
  }

  template<typename T>
  static std::enable_if_t<
    std::is_unsigned<T>::value,
    T
  > read(std::ifstream& file) {
    std::bitset<8 * sizeof(T)> bits;
    file >> bits;
    return static_cast<T>(
      bits.to_ulong()
    );
  }

  static bool canRead(const std::string& filename);

  static void write(
    const std::string& filename,
    const BinaryType& binary
  );

  static BinaryType read(const std::string& filename);
};

InterpretResult interpret(const FileHandler::RawData& data);

} // namespace FileHandlers

} // namespace IO

} // namespace molassembler

#endif
