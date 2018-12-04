// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Molecule/AtomEnvironmentHash.h"

#include "boost/range/iterator_range_core.hpp"
#include "boost/range/join.hpp"
#include "chemical_symmetries/Properties.h"

#include "molassembler/StereopermutatorList.h"
#include "molassembler/Graph/InnerGraph.h"

#include "temple/Functional.h"

namespace molassembler {

namespace hashes {

BondInformation::BondInformation(
  BondType passBondType,
  bool passStereopermutatorOnBond,
  boost::optional<unsigned> passAssignmentOptional
) : bondType(passBondType),
    stereopermutatorOnBond(passStereopermutatorOnBond),
    assignmentOptional(std::move(passAssignmentOptional))
{}

WideHashType BondInformation::hash() const {
  // Initialize with the bond type
  WideHashType hash (
    /* BondType is unsigned, and we want it on the left of the assignment
     * optional, so we shift it by three
     */
    (
      /* The bond type underlying value has to be incremented because otherwise
       * Single is the same as None
       */
      static_cast<
        std::underlying_type_t<BondType>
      >(bondType) + 1
    ) << 3
  );

  /* In the other three bits of given width, we have to store the following
   * cases:
   *
   * - 0: No BondStereopermutator
   * - 1: Unassigned BondStereopermutator
   * - 2-7: BondStereopermutator assignment
   *
   * We can store 6 BondStereopermutator assignments, which ought to be okay.
   * The most you can probably get right now is 5 by fusing something
   * pentagonal at an axial position.
   */
  if(stereopermutatorOnBond) {
    if(assignmentOptional == boost::none) {
      // Information that there is an unassigned BondStereopermutator
      hash += 1;
    } else {
      // Explicit assignment information
      hash += 2 + *assignmentOptional;
    }
  }

  // The remaining case, no stereopermutator on bond, is just hash += 0

  return hash;
}

bool BondInformation::operator < (const BondInformation& other) const {
  return (
    std::tie(bondType, stereopermutatorOnBond, assignmentOptional)
    < std::tie(other.bondType, other.stereopermutatorOnBond, other.assignmentOptional)
  );
}

bool BondInformation::operator == (const BondInformation& other) const {
  return (
    std::tie(bondType, stereopermutatorOnBond, assignmentOptional)
    == std::tie(other.bondType, other.stereopermutatorOnBond, other.assignmentOptional)
  );
}

WideHashType atomEnvironment(
  const temple::Bitmask<AtomEnvironmentComponents>& bitmask,
  const Scine::Utils::ElementType elementType,
  const std::vector<BondInformation>& sortedBonds,
  boost::optional<Symmetry::Name> symmetryNameOptional,
  boost::optional<unsigned> assignedOptional
) {
  static_assert(
    // Sum of information in bits we want to pack in the 128 bit hash
    (
      // Element type (fixed as this cannot possibly increase)
      7
      // Bond information: exactly as many as the largest possible symmetry
      + Symmetry::constexprProperties::maxSymmetrySize * (
        BondInformation::hashWidth
      )
      // The bits needed to store the symmetry name
      + temple::Math::ceil(
        temple::Math::log(Symmetry::nSymmetries + 1.0, 2.0)
      )
      // Roughly 5040 possible assignment values (maximally asymmetric square antiprismatic)
      + 13
    ) <= 128,
    "Element type, bond and symmetry information no longer fit into a 64-bit unsigned integer"
  );

  /* First 7 bits of the 64 bit unsigned number are from the element type
   *
   * Biggest is Cn, which has a value of 112 -> Fits in 7 bits (2^7 = 128)
   */
  WideHashType value;
  if(bitmask & AtomEnvironmentComponents::ElementTypes) {
    value = static_cast<WideHashType>(elementType);
  } else {
    value = 0;
  }

  /* Bond types have 7 possible values currently, plus None is 8
   * -> fits into 3 bits (2^3 = 8).
   *
   * I think bond stereopermutator assignments can be up to 5 (if fused axially
   * onto some pentagonal structure) maximally, but we can fit up to 7 into 3
   * bits.
   *
   * So, a single bond's information needs 6 bits.
   *
   * So, left shift by 7 bits (so there are 7 zeros on the right in the bit
   * representation, where the element type is stored) plus the current bond
   * number multiplied by the width of a BondInformation hash to place a
   * maximum of 8 bond types (maximum symmetry size currently)
   *
   * This occupies 6 * 8 = 48 bits.
   */
  if(bitmask & AtomEnvironmentComponents::BondOrders) {
    unsigned bondNumber = 0;
    for(const auto& bond : sortedBonds) {
      value += (
        bond.hash()
      ) << (7 + BondInformation::hashWidth * bondNumber);

      ++bondNumber;
    }
  }

  if((bitmask & AtomEnvironmentComponents::Symmetries) && symmetryNameOptional) {
    /* We add symmetry information on non-terminal atoms. There are currently
     * 16 symmetries, plus None is 17, which fits into 5 bits (2^5 = 32)
     */
    value += (WideHashType(symmetryNameOptional.value()) + 1) << (7 + 48);

    if(bitmask & AtomEnvironmentComponents::Stereopermutations) {
      /* The remaining space 64 - (7 + 32 + 5) = 20 bits is used for the current
       * permutation. In that space, we can store up to 2^20 - 2 > 1e5
       * permutations (one (0) for no stereopermutator, one (1) for unassigned),
       * which ought to be plenty of bit space. Maximally asymmetric square
       * antiprismatic has around 6k permutations, which fits into 13 bits.
       */
      WideHashType permutationValue;
      if(assignedOptional) {
        permutationValue = assignedOptional.value() + 2;
      } else {
        permutationValue = 1;
      }
      value += permutationValue << (7 + 48 + 5);
    }
  }

  return value;
}

std::vector<WideHashType> generate(
  const InnerGraph& inner,
  const StereopermutatorList& stereopermutators,
  const temple::Bitmask<AtomEnvironmentComponents>& bitmask
) {
  std::vector<WideHashType> hashes;

  const AtomIndex N = inner.N();
  hashes.reserve(N);

  std::vector<BondInformation> bonds;
  bonds.reserve(Symmetry::constexprProperties::maxSymmetrySize);

  for(AtomIndex i = 0; i < N; ++i) {
    if(bitmask & AtomEnvironmentComponents::BondOrders) {
      bonds.clear();

      for(
        const InnerGraph::Edge& edge :
        boost::make_iterator_range(inner.edges(i))
      ) {
        auto stereopermutatorOption = stereopermutators.option(
          BondIndex {
            inner.source(edge),
            inner.target(edge)
          }
        );

        if(stereopermutatorOption) {
          bonds.emplace_back(
            inner.bondType(edge),
            true,
            stereopermutatorOption->assigned()
          );
        } else {
          bonds.emplace_back(
            inner.bondType(edge),
            false,
            boost::none
          );
        }
      }

      std::sort(
        bonds.begin(),
        bonds.end()
      );
    }

    boost::optional<Symmetry::Name> symmetryNameOption;
    boost::optional<unsigned> assignmentOption;

    if(auto refOption = stereopermutators.option(i)) {
      symmetryNameOption = refOption->getSymmetry();
      assignmentOption = refOption->assigned();
    }

    hashes.emplace_back(
      atomEnvironment(
        bitmask,
        inner.elementType(i),
        bonds,
        symmetryNameOption,
        assignmentOption
      )
    );
  }

  return hashes;
}

std::tuple<
  std::vector<HashType>,
  std::vector<HashType>,
  HashType
> narrow(
  const std::vector<WideHashType>& a,
  const std::vector<WideHashType>& b
) {
  std::unordered_map<WideHashType, HashType> reductionMapping;
  HashType counter = 0;

  // Generate mapping from hash values to integer-incremented reduction
  for(const auto& hash : boost::range::join(a, b)) {
    if(reductionMapping.count(hash) == 0) {
      reductionMapping.emplace(
        hash,
        counter
      );

      ++counter;
    }
  }

  // counter now contains the maximum new hash for this set of hashes
  return {
    temple::map(
      a,
      [&reductionMapping](const WideHashType& hash) -> HashType {
        return reductionMapping.at(hash);
      }
    ),
    temple::map(
      b,
      [&reductionMapping](const WideHashType& hash) -> HashType {
        return reductionMapping.at(hash);
      }
    ),
    counter
  };
}

/* Ensure that the LookupFunctor matches the required boost concept. Although
 * the documentation asks for the UnaryFunction concept, in reality it requires
 * AdaptableUnaryFunction with typedefs for argument and result.
 */
BOOST_CONCEPT_ASSERT((
  boost::AdaptableUnaryFunctionConcept<LookupFunctor, WideHashType, AtomIndex>
));

} // namespace hashes

} // namespace molassembler
