#ifndef INCLUDE_PERSISTENT_INDEX_MOLECULE_INTERFACE_H
#define INCLUDE_PERSISTENT_INDEX_MOLECULE_INTERFACE_H

#include "IO.h"

/*! @file
 *
 * Defines a proxy class that aids in atom index book-keeping. 
 */

namespace molassembler {

/*!
 * Book-keeping helper class that keeps track of deleted atom indices so that
 * all indices are persistent.
 */
class PersistentIndicesInterface {
private:
  std::unique_ptr<Molecule> _molecule;
  std::set<AtomIndexType> _deletedIndices;
  /* Converts an external persistent index into the minimal internal one. This
   * throws if either the internal index has been deleted or the resulting 
   * internal index is out of bounds (i.e. larger than the amount of internal
   * atoms).
   */
  AtomIndexType _toInternalIndex(const AtomIndexType& a) const;

  /* Converts an internal index into the persistent external index. Throws if 
   * the internal index passed to it is out of bounds (i.e. larger than the 
   * biggest internal index).
   */
  AtomIndexType _toExternalIndex(const AtomIndexType& a) const;

public:
  // Constructors from file, two elements and shared bond
  PersistentIndicesInterface();
  PersistentIndicesInterface(const std::string& molFile);
  PersistentIndicesInterface(
    const Delib::ElementType& a,
    const Delib::ElementType& b,
    const BondType& bondType
  );

  // Mirror entire public modification interface to Molecule
  // Maybe with templates? Hana?
  // Rewriting every single function with forwarding calls isn't so attractive
};

PersistentIndicesInterface::PersistentIndicesInterface() {} 
PersistentIndicesInterface::PersistentIndicesInterface(
  const std::string& molFile
) {
  IO::MOLFileHandler molReader;
  _molecule = std::make_unique<Molecule>(
    molReader.readSingle(molFile)
  );
} 

PersistentIndicesInterface::PersistentIndicesInterface(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) {
  _molecule = std::make_unique<Molecule>(
    a,
    b,
    bondType
  );
}

AtomIndexType PersistentIndicesInterface::_toInternalIndex(const AtomIndexType& a) const {
  AtomIndexType returnIndex = a;

  for(const auto& deletedExternalIndex : _deletedIndices) {
    if(deletedExternalIndex > a) returnIndex--;
  }

  if(
    _deletedIndices.count(a) != 0
    || returnIndex < _molecule->getNumAtoms()
  ) {
    throw std::logic_error(
      "Requested internal index is either deleted or out of bounds."
    );
  }

  return returnIndex;
}

AtomIndexType PersistentIndicesInterface::_toExternalIndex(const AtomIndexType& a) const {
  if(a >= _molecule->getNumAtoms()) {
    throw std::logic_error(
      "Conversion to external index requested for out of bounds internal index."
    );
  }

  AtomIndexType returnIndex = a;

  for(const auto& deletedExternalIndex : _deletedIndices) {
    if(deletedExternalIndex <= a) returnIndex++;
  }

  return returnIndex;
}

} // namespace molassembler

#endif

