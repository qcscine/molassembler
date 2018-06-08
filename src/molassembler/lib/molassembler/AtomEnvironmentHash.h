#ifndef INCLUDE_MOLASSEMBLER_ATOM_ENVIRONMENT_HASH_H
#define INCLUDE_MOLASSEMBLER_ATOM_ENVIRONMENT_HASH_H

#include "chemical_symmetries/Symmetries.h"
#include "temple/constexpr/Bitmask.h"

#include "detail/SharedTypes.h"

/*!@file
 *
 * Contains a function for hashing an atom's environment in a Molecule
 */

namespace molassembler {

namespace hashes {

using AtomEnvironmentHashType = unsigned long long;

//! Convolutes the atom's element type and bonds into a characteristic number
AtomEnvironmentHashType atomEnvironment(
  const temple::Bitmask<AtomEnvironmentComponents>& bitmask,
  const Delib::ElementType elementType,
  const std::vector<BondType>& sortedBonds,
  boost::optional<Symmetry::Name> symmetryNameOptional,
  boost::optional<unsigned> assignedOptional
);

} // namespace hashes

} // namespace molassembler

#endif
