#ifndef INCLUDE_MOLASSEMBLER_ATOM_ENVIRONMENT_HASH_H
#define INCLUDE_MOLASSEMBLER_ATOM_ENVIRONMENT_HASH_H

#include "boost/multiprecision/cpp_int.hpp"
#include "boost/optional.hpp"
#include "Delib/ElementTypes.h"

#include "chemical_symmetries/Names.h"
#include "temple/constexpr/Bitmask.h"

#include "molassembler/Types.h"

#include <vector>

/*!@file
 *
 * Contains a function for hashing an atom's environment in a Molecule
 */

namespace molassembler {

// Forward-declarations
class StereocenterList;
class InnerGraph;

namespace hashes {

using WideHashType = boost::multiprecision::uint128_t;
using HashType = std::uint64_t;

struct BondInformation {
  static constexpr unsigned hashWidth = 6;

  BondType bondType;
  bool stereocenterOnBond;
  boost::optional<unsigned> assignmentOptional;

  BondInformation(
    BondType passBondType,
    bool passStereocenterOnBond,
    boost::optional<unsigned> passAssignmentOptional
  );

  WideHashType hash() const;

  bool operator < (const BondInformation& other) const;
  bool operator == (const BondInformation& other) const;
};

//! Convolutes the atom's element type and bonds into a characteristic number
WideHashType atomEnvironment(
  const temple::Bitmask<AtomEnvironmentComponents>& bitmask,
  Delib::ElementType elementType,
  const std::vector<BondInformation>& sortedBonds,
  boost::optional<Symmetry::Name> symmetryNameOptional,
  boost::optional<unsigned> assignedOptional
);

//! Generates the hashes for every atom in a molecule's components
std::vector<WideHashType> generate(
  const InnerGraph& inner,
  const StereocenterList& stereocenters,
  const temple::Bitmask<AtomEnvironmentComponents>& bitmask
);

//! Re-enumerates the hashes in two generated hash lists
std::tuple<
  std::vector<HashType>,
  std::vector<HashType>,
  HashType
> narrow(
  const std::vector<WideHashType>& a,
  const std::vector<WideHashType>& b
);

std::pair<
  std::vector<HashType>,
  std::vector<HashType>
> generate(
  const InnerGraph& aGraph,
  const StereocenterList& aStereocenters,
  const InnerGraph& bGraph,
  const StereocenterList& bStereocenters,
  const temple::Bitmask<AtomEnvironmentComponents>& bitmask
);

/*! A functor for getting an atom's hash from a captured list of hashes
 *
 * \note This explicit form is needed instead of a lambda for boost's concept
 *   checks. It satisfies the AdaptableUnaryFunctionConcept.
 */
struct LookupFunctor {
  const std::vector<HashType>* const hashesPtr;

  using argument_type = AtomIndex;
  using result_type = HashType;

  inline LookupFunctor(const std::vector<HashType>& hashes) : hashesPtr(&hashes) {}

  inline HashType operator() (const AtomIndex i) const {
    return hashesPtr->at(i);
  }
};

} // namespace hashes

} // namespace molassembler

#endif
