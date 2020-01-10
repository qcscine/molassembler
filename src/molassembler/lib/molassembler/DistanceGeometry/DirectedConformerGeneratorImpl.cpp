/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/DistanceGeometry/DirectedConformerGeneratorImpl.h"

#include "molassembler/Cycles.h"
#include "molassembler/Molecule/MoleculeImpl.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/StereopermutatorList.h"

#include "stereopermutation/Composites.h"
#include "temple/Functional.h"

#include "Utils/Geometry/AtomCollection.h"

#include "boost/variant.hpp"

namespace Scine {
namespace molassembler {

namespace detail {

//! C++'s modulo yields sign of the dividend, but we want nonnegative-always
int euclideanModulo(const int a, const int base) {
  return ((a % base) + base) % base;
}

//! Distance between elements in a modulus value cycle (0 -> 1 -> 2 -> 0)
unsigned distance(const int i, const int j, const int U) {
  assert(i >= 0 && j >= 0 && U > 1);

  int value = std::min(
    euclideanModulo(i - j, U),
    euclideanModulo(j - i, U)
  );

  assert(value > 0);

  return static_cast<unsigned>(value);
}

template<typename ChoiceIndex>
struct BoundedNodeTrieChooseFunctor {
  /*!
   * @brief Calculates merit of a particular choice in a modulus value cycle
   *   with some elements in the cycle existing and others not.
   */
  double merit(const ChoiceIndex choice, const boost::dynamic_bitset<>& childExists) {
    double value = 0.0;

    const ChoiceIndex U = childExists.size();

    for(ChoiceIndex other = 0; other < U; ++other) {
      if(other == choice) {
        continue;
      }

      if(childExists.test(other)) {
        unsigned choiceDistance = distance(choice, other, U);
        assert(choiceDistance < std::numeric_limits<ChoiceIndex>::max());
        value += choiceDistance;
      }
    }

    assert(value >= 0);

    return value;
  }

  /*!
   * @brief At every depth, try to make a locally optimal choice
   *
   * Of the viable children, pick that one which has highest merit, i.e. is
   * most distant from all already existing children.
   */
  ChoiceIndex operator () (const std::vector<ChoiceIndex>& viableChildren, const boost::dynamic_bitset<>& existingChildren) {
    assert(!viableChildren.empty());

    /* In case not all children already exist, it makes sense to calculate
     * merits for all choices at this level.
     */
    if(!existingChildren.all()) {
      std::vector<ChoiceIndex> bestChoices;
      double bestMerit = 0;
      for(ChoiceIndex viableChild : viableChildren) {
        if(!existingChildren.test(viableChild)) {
          double choiceMerit = merit(viableChild, existingChildren);
          if(choiceMerit > bestMerit) {
            bestChoices = {viableChild};
            bestMerit = choiceMerit;
          } else if(choiceMerit == bestMerit) {
            bestChoices.push_back(viableChild);
          }
        }
      }

      assert(!bestChoices.empty());

      return bestChoices.at(
        std::uniform_int_distribution<ChoiceIndex>(
          0,
          bestChoices.size() - 1
        )(randomnessEngine())
      );
    }

    /* All children exist, so no point in calculating merits: we choose from
     * the viable ones
     */
    return viableChildren.at(
      std::uniform_int_distribution<ChoiceIndex>(0, viableChildren.size() - 1)(randomnessEngine())
    );
  }
};

} // namespace detail

unsigned DirectedConformerGenerator::Impl::distance(
  const DecisionList& a,
  const DecisionList& b,
  const DecisionList& bounds
) {
  if(a.size() != b.size() || a.size() != bounds.size()) {
    throw std::invalid_argument("Not all decision lists have the same length");
  }

  unsigned distance = 0;
  for(unsigned i = 0; i < bounds.size(); ++i) {
    assert(a[i] < bounds[i] && b[i] < bounds[i]);
    distance += detail::distance(a[i], b[i], bounds[i]);
  }

  return distance;
}


boost::variant<DirectedConformerGenerator::IgnoreReason, BondStereopermutator>
DirectedConformerGenerator::Impl::considerBond(
  const BondIndex& bondIndex,
  const Molecule& molecule,
  const std::unordered_map<AtomIndex, unsigned>& smallestCycleMap
) {
  // Make sure the bond exists in the first place
  assert(molecule.graph().adjacent(bondIndex.first, bondIndex.second));

  // If either atom is terminal, we can stop immediately
  if(
    molecule.graph().degree(bondIndex.first) == 1
    || molecule.graph().degree(bondIndex.second) == 1
  ) {
    return IgnoreReason::HasTerminalConstitutingAtom;
  }

  BondType bondType = molecule.graph().bondType(bondIndex);
  if(bondType == BondType::Eta) {
    return IgnoreReason::IsEtaBond;
  }

  /* The dihedral already has an assigned bond stereopermutator placed on it or
   * that bond stereopermutator is isotropic (ranking of all substituents at
   * one side is the same)
   */
  if(auto bondStereopermutatorOption = molecule.stereopermutators().option(bondIndex)) {
    if(bondStereopermutatorOption->assigned()) {
      return IgnoreReason::HasAssignedBondStereopermutator;
    }

    if(bondStereopermutatorOption->composite().isIsotropic()) {
      return IgnoreReason::RotationIsIsotropic;
    }
  }

  // Check with cycles information next: If the dihedral is in a cycle, skip it
  if(
    smallestCycleMap.count(bondIndex.first) > 0
    && smallestCycleMap.count(bondIndex.second) > 0
  ) {
    return IgnoreReason::InCycle;
  }

  /* Last, instantiate a BondStereopermutator on the bond and use that to
   * figure out rotational isotropicity
   */
  auto firstAtomStereopermutatorOption = molecule.stereopermutators().option(bondIndex.first);
  auto secondAtomStereopermutatorOption = molecule.stereopermutators().option(bondIndex.second);

  // Both atom stereopermutators the bond will constitute have to exist
  if(!firstAtomStereopermutatorOption || !secondAtomStereopermutatorOption) {
    return IgnoreReason::AtomStereopermutatorPreconditionsUnmet;
  }

  // Try instantiating the BondStereopermutator
  try {
    auto alignment = BondStereopermutator::Alignment::Eclipsed;
    if(bondType == BondType::Single) {
      alignment = BondStereopermutator::Alignment::Staggered;
    }

    BondStereopermutator trialPermutator {
      *firstAtomStereopermutatorOption,
      *secondAtomStereopermutatorOption,
      bondIndex,
      alignment
    };

    if(trialPermutator.numStereopermutations() == 1) {
      return IgnoreReason::RotationIsIsotropic;
    }

    return trialPermutator;
  } catch(...) {
    return IgnoreReason::AtomStereopermutatorPreconditionsUnmet;
  }
}

DirectedConformerGenerator::Impl::Impl(
  Molecule molecule,
  const BondList& bondsToConsider
) : _molecule(std::move(molecule))
{
  // Precalculate the smallest cycle map
  auto smallestCycleMap = makeSmallestCycleMap(_molecule.graph().cycles());

  // Allocate some space for the relevant bonds and new stereopermutators
  _relevantBonds.reserve(_molecule.graph().B() / 2);
  std::vector<BondStereopermutator> bondStereopermutators;
  bondStereopermutators.reserve(_molecule.graph().B() / 2);

  // Make a visitor to move out bond stereopermutators from the variant
  struct ExtractStereopermutatorVisitor : public boost::static_visitor<bool> {
    ExtractStereopermutatorVisitor(std::vector<BondStereopermutator>& stereoList) : vectorRef(stereoList) {}

    std::vector<BondStereopermutator>& vectorRef;

    bool operator() (BondStereopermutator& permutator) {
      vectorRef.push_back(std::move(permutator));
      return true;
    }

    bool operator() (IgnoreReason /* reason */) { return false; }
  };

  ExtractStereopermutatorVisitor visitor {bondStereopermutators};

  if(bondsToConsider.empty()) {
    // Consider all bonds
    for(BondIndex bondIndex : boost::make_iterator_range(_molecule.graph().bonds())) {
      auto importanceVariant = considerBond(bondIndex, _molecule, smallestCycleMap);

      if(boost::apply_visitor(visitor, importanceVariant)) {
        _relevantBonds.push_back(bondIndex);
      }
    }
  } else {
    for(const BondIndex& bondIndex : bondsToConsider) {
      auto importanceVariant = considerBond(bondIndex, _molecule, smallestCycleMap);

      if(boost::apply_visitor(visitor, importanceVariant)) {
        _relevantBonds.push_back(bondIndex);
      }
    }
  }

  // Sort the relevant bonds and shrink
  temple::inplace::sort(_relevantBonds);
  _relevantBonds.shrink_to_fit();

  /* In case there are no bonds to consider, then we're done. No other work
   * to be done.
   */
  if(_relevantBonds.empty()) {
    return;
  }

  // Add the BondStereopermutators to our underlying molecule's stereopermutator list
  for(auto& bondStereopermutator : bondStereopermutators) {
    _molecule._pImpl->_stereopermutators.add(
      std::move(bondStereopermutator)
    );
  }

  // Collect the bounds on each stereopermutator's permutations for the trie
  auto bounds = temple::map(
    _relevantBonds,
    [&](const BondIndex& bond) -> std::uint8_t {
      auto stereoOption = _molecule.stereopermutators().option(bond);
      if(!stereoOption) {
        throw std::logic_error(
          "No BondStereopermutator for relevant bond, but it was just inserted"
        );
      }

      unsigned nPermutations = stereoOption->numStereopermutations();
      assert(nPermutations < std::numeric_limits<std::uint8_t>::max());
      return nPermutations;
    }
  );

  // Initialize the trie with the list of bounds
  _decisionLists.setBounds(std::move(bounds));
}

DirectedConformerGenerator::DecisionList
DirectedConformerGenerator::Impl::generateNewDecisionList() {
  if(_relevantBonds.empty()) {
    throw std::logic_error("List of relevant bonds is empty!");
  }

  detail::BoundedNodeTrieChooseFunctor<std::uint8_t> chooseFunctor {};
  return _decisionLists.generateNewEntry(chooseFunctor);
}

const Molecule& DirectedConformerGenerator::Impl::conformationMolecule(const DecisionList& decisionList) {
  const unsigned U = decisionList.size();

  if(U != _decisionLists.bounds().size() || U != _relevantBonds.size()) {
    throw std::invalid_argument("Passed decision list has wrong length");
  }

  for(unsigned i = 0; i < U; ++i) {
    const BondIndex& bondIndex = _relevantBonds.at(i);
    auto stereoOption = _molecule._pImpl->_stereopermutators.option(bondIndex);
    assert(stereoOption);
    stereoOption->assign(decisionList[i]);
  }

  return _molecule;
}

outcome::result<Utils::PositionCollection>
DirectedConformerGenerator::Impl::generateRandomConformation(
  const DecisionList& decisionList,
  const distance_geometry::Configuration& configuration
) {
  return Scine::molassembler::generateRandomConformation(
    conformationMolecule(decisionList),
    configuration
  );
}

outcome::result<Utils::PositionCollection>
DirectedConformerGenerator::Impl::generateConformation(
  const DecisionList& decisionList,
  const unsigned seed,
  const distance_geometry::Configuration& configuration
) {
  return Scine::molassembler::generateConformation(
    conformationMolecule(decisionList),
    seed,
    configuration
  );
}

DirectedConformerGenerator::DecisionList
DirectedConformerGenerator::Impl::getDecisionList(
  const Utils::AtomCollection& atomCollection,
  const BondStereopermutator::FittingMode mode
) {
  if(
    !temple::all_of(
      boost::make_iterator_range(_molecule.graph().atoms()),
      [&](const AtomIndex i) -> bool {
        return atomCollection.getElement(i) == _molecule.graph().elementType(i);
      }
    )
  ) {
    throw std::logic_error("Input AtomCollection elements do not match generator's underlying molecule. Misordered? Different molecule input?");
  }

  return getDecisionList(atomCollection.getPositions(), mode);
}

DirectedConformerGenerator::DecisionList
DirectedConformerGenerator::Impl::getDecisionList(
  const Utils::PositionCollection& positions,
  const BondStereopermutator::FittingMode mode
) {
  AngstromWrapper angstromPositions {positions};

  if(
    !temple::all_of(
      _relevantBonds,
      [&](const BondIndex& bondIndex) -> bool {
        const auto& permutators = _molecule._pImpl->_stereopermutators;
        return (
          permutators.option(bondIndex.first)
          && permutators.option(bondIndex.second)
          && permutators.option(bondIndex)
        );
      }
    )
  ) {
    throw std::logic_error("Underlying molecule permutator preconditions unmet!");
  }

  return temple::map(
    _relevantBonds,
    [&](const BondIndex& bondIndex) -> std::uint8_t {
      auto firstAtom = _molecule._pImpl->_stereopermutators.option(bondIndex.first);
      auto secondAtom = _molecule._pImpl->_stereopermutators.option(bondIndex.second);
      auto stereoOption = _molecule._pImpl->_stereopermutators.option(bondIndex);

      assert(firstAtom && secondAtom && stereoOption);
      stereoOption->fit(angstromPositions, *firstAtom, *secondAtom, mode);

      return stereoOption->assigned().value_or(unknownDecision);
    }
  );
}

} // namespace molassembler
} // namespace Scine
