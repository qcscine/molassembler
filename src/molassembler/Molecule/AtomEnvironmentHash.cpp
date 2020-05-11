/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Molecule/AtomEnvironmentHash.h"

#include "boost/range/join.hpp"
#include "molassembler/Shapes/PropertyCaching.h"

#include "molassembler/AtomStereopermutator.h"
#include "molassembler/BondStereopermutator.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/Graph/PrivateGraph.h"

#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/Adaptors/Iota.h"

#include "Utils/Geometry/ElementInfo.h"

namespace Scine {
namespace Molassembler {
namespace Hashes {

BondInformation::BondInformation(
  BondType passBondType,
  bool passStereopermutatorOnBond,
  boost::optional<unsigned> passAssignmentOptional
) : bondType(passBondType),
    stereopermutatorOnBond(passStereopermutatorOnBond),
    assignmentOptional(std::move(passAssignmentOptional))
{}

static_assert(
  BondInformation::bondTypeBits == Temple::Math::ceil(
    Temple::Math::log(nBondTypes + 1.0, 2.0)
  ),
  "Number of bond types requires a change in bond type bits!"
);

static_assert(
  BondInformation::hashWidth == BondInformation::bondTypeBits + BondInformation::assignmentBits,
  "BondInformation::hashWidth is no longer the sum of bond type bits and assignment bits"
);

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
    ) << bondTypeBits
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

  // The remaining case, no stereopermutator on bond, is just hash = 0
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

WideHashType hash(
  const AtomEnvironmentComponents bitmask,
  const Utils::ElementType elementType,
  const std::vector<BondInformation>& sortedBonds,
  const boost::optional<Shapes::Shape>& shapeOptional,
  const boost::optional<unsigned>& assignedOptional
) {
  /* The bit representation of element types is 16 bits wide storing both atomic
   * number and atomic mass number in order to account for isotopes
   */
  constexpr unsigned elementTypeBits = 16;
  constexpr unsigned shapeNameBits = Temple::Math::ceil(
    Temple::Math::log(Shapes::nShapes + 1.0, 2.0)
  );

  static_assert(
    // Sum of information in bits we want to pack in the 128 bit hash
    (
      // Element type (fixed as this cannot possibly increase)
      elementTypeBits
      // Bond information: exactly as many as the largest possible shape
      + Shapes::ConstexprProperties::maxShapeSize * (
        BondInformation::hashWidth
      )
      // The bits needed to store the shape name (plus none)
      + shapeNameBits
      // Roughly 5040 possible assignment values (maximally asymmetric square antiprismatic)
      + 13
    ) <= 128,
    "Element type, bond and shape information no longer fit into a 64-bit unsigned integer"
  );

  // First 16 bits of the number are from the element type
  WideHashType value;
  if(bitmask & AtomEnvironmentComponents::ElementTypes) {
    using ElementTypeUnderlying = std::underlying_type<Utils::ElementType>::type;
    value = static_cast<WideHashType>(static_cast<ElementTypeUnderlying>(elementType));
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
   * So, left shift by 16 bits (so there are 16 zeros on the right in the bit
   * representation, where the element type is stored) plus the current bond
   * number multiplied by the width of a BondInformation hash to place a
   * maximum of 12 bond types (maximum shape size currently)
   *
   * This occupies 6 * 12 = 72 bits.
   */
  constexpr unsigned bondsHashSectionWidth = BondInformation::hashWidth * Shapes::ConstexprProperties::maxShapeSize;

  if(bitmask & AtomEnvironmentComponents::BondOrders) {
    unsigned bondNumber = 0;
    for(const auto& bond : sortedBonds) {
      value += (
        bond.hash()
      ) << (elementTypeBits + BondInformation::hashWidth * bondNumber);

      ++bondNumber;
    }
  }

  if((bitmask & AtomEnvironmentComponents::Shapes) && shapeOptional) {
    /* We add shape information on non-terminal atoms. There are currently
     * 30 shapes, plus None is 31, which fits into 5 bits (2^5 = 32)
     */
    value += (WideHashType(shapeOptional.value()) + 1) << (elementTypeBits + bondsHashSectionWidth);

    if(bitmask & AtomEnvironmentComponents::Stereopermutations) {
      /* The remaining space (128 - (16 + 72 + 5) = 35 bits) is used for the
       * current permutation. Log_2(12!) ~= 29, so we're good.
       */
      WideHashType permutationValue;
      if(assignedOptional) {
        permutationValue = assignedOptional.value() + 2;
      } else {
        permutationValue = 1;
      }
      value += permutationValue << (elementTypeBits + bondsHashSectionWidth + shapeNameBits);
    }
  }

  return value;
}

std::vector<BondInformation> gatherBonds(
  const PrivateGraph& inner,
  boost::optional<const StereopermutatorList&> stereopermutators,
  const AtomEnvironmentComponents componentsBitmask,
  const AtomIndex i
) {
  std::vector<BondInformation> bonds;
  bonds.reserve(Shapes::ConstexprProperties::maxShapeSize);

  if(componentsBitmask & AtomEnvironmentComponents::Stereopermutations) {
    for(const PrivateGraph::Edge& edge : inner.edges(i)) {
      const BondIndex bond {
        inner.source(edge),
        inner.target(edge)
      };

      if(stereopermutators) {
        auto stereopermutatorOption = stereopermutators->option(bond);

        if(stereopermutatorOption && stereopermutatorOption->numAssignments() > 1) {
          bonds.emplace_back(
            inner.bondType(edge),
            true,
            stereopermutatorOption->assigned()
          );
        } else {
          /* Even if a stereopermutator is present, if it has only a single
           * viable assignment, it is best that it cannot contribute to
           * differentiating between molecules. This avoids e.g. that spuriously
           * different interpretations of aromatic cycles due to
           * close-to-threshold bond ordes lead to false distinctions between
           * molecules.
           */
          bonds.emplace_back(
            inner.bondType(edge),
            false,
            boost::none
          );
        }
      }
    }
  } else {
    for(const PrivateGraph::Edge& edge : inner.edges(i)) {
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

  return bonds;
}

WideHashType atomEnvironment(
  const PrivateGraph& inner,
  boost::optional<const StereopermutatorList&> stereopermutators,
  AtomEnvironmentComponents bitmask,
  AtomIndex i
) {
  std::vector<BondInformation> bonds;
  boost::optional<Shapes::Shape> shapeOption;
  boost::optional<unsigned> assignmentOption;

  if(bitmask & AtomEnvironmentComponents::BondOrders) {
    bonds = gatherBonds(inner, stereopermutators, bitmask, i);
  }

  if(stereopermutators) {
    if(auto refOption = stereopermutators->option(i)) {
      shapeOption = refOption->getShape();
      assignmentOption = refOption->assigned();
    }
  }

  return hash(
    bitmask,
    inner.elementType(i),
    bonds,
    shapeOption,
    assignmentOption
  );
}

std::vector<WideHashType> generate(
  const PrivateGraph& inner,
  boost::optional<const StereopermutatorList&> stereopermutators,
  const AtomEnvironmentComponents bitmask
) {
  const unsigned N = inner.N();
  std::vector<WideHashType> hashes(N);

#pragma omp parallel for
  for(unsigned i = 0; i < N; ++i) {
    hashes.at(i) = atomEnvironment(
      inner,
      stereopermutators,
      bitmask,
      i
    );
  }

  return hashes;
}

bool identityCompare(
  const PrivateGraph& aGraph,
  const StereopermutatorList& aStereopermutators,
  const PrivateGraph& bGraph,
  const StereopermutatorList& bStereopermutators,
  AtomEnvironmentComponents componentBitmask
) {
  assert(aGraph.N() == bGraph.N());

  return Temple::all_of(
    Temple::Adaptors::range(aGraph.N()),
    [&](const AtomIndex i) -> WideHashType {
      return atomEnvironment(
        aGraph,
        aStereopermutators,
        componentBitmask,
        i
      ) == atomEnvironment(
        bGraph,
        bStereopermutators,
        componentBitmask,
        i
      );
    }
  );
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
    Temple::map(
      a,
      [&reductionMapping](const WideHashType& hash) -> HashType {
        return reductionMapping.at(hash);
      }
    ),
    Temple::map(
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

} // namespace Hashes
} // namespace Molassembler
} // namespace Scine
