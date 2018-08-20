#ifndef INCLUDE_MOLASSEMBLER_ATOM_ENVIRONMENT_HASH_H
#define INCLUDE_MOLASSEMBLER_ATOM_ENVIRONMENT_HASH_H

#include "chemical_symmetries/Names.h"
#include "temple/constexpr/Bitmask.h"

#include "molassembler/detail/SharedTypes.h"

/*!@file
 *
 * Contains a function for hashing an atom's environment in a Molecule
 */

namespace molassembler {

// Forward-declarations
class StereocenterList;

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

//! Generates the hashes for every atom in a molecule's components
std::vector<AtomEnvironmentHashType> generate(
  const GraphType& graph,
  const StereocenterList& stereocenters,
  const temple::Bitmask<AtomEnvironmentComponents>& bitmask
);

//! Re-enumerates the hashes in two generated hash lists
AtomEnvironmentHashType regularize(
  std::vector<AtomEnvironmentHashType>& a,
  std::vector<AtomEnvironmentHashType>& b
);

/*! A functor for getting an atom's hash from a captured list of hashes
 *
 * \note This explicit form is needed instead of a lambda for boost's concept
 *   checks. It satisfies the AdaptableUnaryFunctionConcept.
 */
struct LookupFunctor {
  const std::vector<AtomEnvironmentHashType>* const hashes;

  using argument_type = AtomIndexType;
  using result_type = AtomEnvironmentHashType;

  inline LookupFunctor(const std::vector<AtomEnvironmentHashType>& hashes) : hashes(&hashes) {}

  inline AtomEnvironmentHashType operator() (const AtomIndexType i) const {
    return hashes->at(i);
  }
};

} // namespace hashes

} // namespace molassembler

#endif
