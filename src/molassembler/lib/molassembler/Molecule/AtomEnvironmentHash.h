/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Hash an atom's environment in a Molecule for isomorphism calculations
 */

#ifndef INCLUDE_MOLASSEMBLER_ATOM_ENVIRONMENT_HASH_H
#define INCLUDE_MOLASSEMBLER_ATOM_ENVIRONMENT_HASH_H

#include "boost/multiprecision/cpp_int.hpp"
#include "boost/optional.hpp"
#include "Utils/Geometry/ElementTypes.h"

#include "shapes/Shapes.h"

#include "molassembler/Types.h"

#include <vector>

namespace Scine {

namespace molassembler {

// Forward-declarations
class StereopermutatorList;
class InnerGraph;

/**
 * @brief Classes and methods to compute hashes of atom environments
 */
namespace hashes {

using WideHashType = boost::multiprecision::uint128_t;
using HashType = std::uint64_t;

/**
 * @brief Information pertinent to a singular bond hash
 */
struct BondInformation {
  static constexpr unsigned bondTypeBits = 3;
  static constexpr unsigned assignmentBits = 3;
  /* - 3 bits for bondtype + none
   * - 3 bits for the assignment (none, unassigned, values)
   */
  static constexpr unsigned hashWidth = 6;

  BondType bondType;
  bool stereopermutatorOnBond;
  boost::optional<unsigned> assignmentOptional;

  BondInformation(
    BondType passBondType,
    bool passStereopermutatorOnBond,
    boost::optional<unsigned> passAssignmentOptional
  );

  WideHashType hash() const;

  bool operator < (const BondInformation& other) const;
  bool operator == (const BondInformation& other) const;
};

/*! @brief Convolutes the atom's element type and bonds into an unsigned integer
 *
 * At current, this is a bijective mapping and has zero probability of
 * collisions.
 *
 * @complexity{@math{\Theta(1)}}
 *
 * @todo Whether a stereopermutator is a stereocenter or not currently is
 *   not distinguished in the hash (i.e. 0/1 is treated the same as 0/2).
 */
WideHashType hash(
  AtomEnvironmentComponents bitmask,
  Utils::ElementType elementType,
  const std::vector<BondInformation>& sortedBonds,
  const boost::optional<Shapes::Shape>& shapeOptional,
  const boost::optional<unsigned>& assignedOptional
);

/*! @brief Collects bond information for hash() input
 *
 * @complexity{@math{\Theta(1)} since the out-degree of atoms is considered a
 * fixed small number because molecules are sparse.}
 */
std::vector<BondInformation> gatherBonds(
  const InnerGraph& inner,
  const StereopermutatorList& stereopermutators,
  AtomEnvironmentComponents componentsBitmask,
  AtomIndex i
);

/*! @brief Calculate the hash for a particular atom index
 *
 * @complexity{@math{\Theta(1)}}
 */
WideHashType atomEnvironment(
  const InnerGraph& inner,
  const StereopermutatorList& stereopermutators,
  AtomEnvironmentComponents bitmask,
  AtomIndex i
);

/*! @brief Generates the hashes for every atom in a molecule's components
 *
 * @complexity{@math{\Theta(N)}}
 */
std::vector<WideHashType> generate(
  const InnerGraph& inner,
  const StereopermutatorList& stereopermutators,
  AtomEnvironmentComponents bitmask
);

/*! @brief Compare two molecules with an identity index mapping
 *
 * @complexity{@math{O(N)}}
 */
bool identityCompare(
  const InnerGraph& aGraph,
  const StereopermutatorList& aStereopermutators,
  const InnerGraph& bGraph,
  const StereopermutatorList& bStereopermutators,
  AtomEnvironmentComponents componentBitmask
);

/*! @brief Re-enumerates the hashes in two generated hash lists
 *
 * Although the hash space requires a wide integer type, the number of different
 * atom environments in molecules will not require a wide integer type. We can
 * therefore remap the wide hashes to an enumerated smaller integer type.
 *
 * @complexity{@math{\Theta(N)}}
 */
std::tuple<
  std::vector<HashType>,
  std::vector<HashType>,
  HashType
> narrow(
  const std::vector<WideHashType>& a,
  const std::vector<WideHashType>& b
);

/*! @brief Generates hashes for two molecules' components
 *
 * Generate hashes for two molecules. Basis for comparing vertices in an
 * isomorphism algorithm.
 *
 * @complexity{@math{\Theta(N)}}
 */
std::pair<
  std::vector<HashType>,
  std::vector<HashType>
> generate(
  const InnerGraph& aGraph,
  const StereopermutatorList& aStereopermutators,
  const InnerGraph& bGraph,
  const StereopermutatorList& bStereopermutators,
  AtomEnvironmentComponents bitmask
);

/*! @brief A functor for getting an atom's hash from a captured list of hashes
 *
 * @note This explicit form is needed instead of a lambda for boost's concept
 *   checks. It must satisfy AdaptableUnaryFunctionConcept.
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

} // namespace Scine

#endif
