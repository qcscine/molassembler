/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "Molassembler/DistanceGeometry/DirectedConformerGeneratorImpl.h"

#include "Molassembler/Cycles.h"
#include "Molassembler/Molecule/MoleculeImpl.h"
#include "Molassembler/Graph.h"
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/Detail/Cartesian.h"

#include "Molassembler/Stereopermutation/Composites.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Adaptors/SequentialPairs.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Optionals.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/DistanceGeometry/Error.h"

#include "Utils/Geometry/AtomCollection.h"
#include "boost/variant.hpp"

namespace Scine {
namespace Molassembler {
namespace Detail {

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
  Random::Engine& engine;

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
      return Temple::Random::pick(bestChoices, engine);
    }

    // All children exist, so no point in calculating merits
    return Temple::Random::pick(viableChildren, engine);
  }
};

} // namespace Detail

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
    distance += Detail::distance(a[i], b[i], bounds[i]);
  }

  return distance;
}


boost::variant<DirectedConformerGenerator::IgnoreReason, BondStereopermutator>
DirectedConformerGenerator::Impl::considerBond(
  const BondIndex& bondIndex,
  const Molecule& molecule,
  BondStereopermutator::Alignment alignment
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

  // Check with cycles information next: If the bond is in a cycle, skip it
  if(molecule.graph().cycles().numCycleFamilies(bondIndex) > 0) {
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
    if(
      alignment != BondStereopermutator::Alignment::EclipsedAndStaggered
      && alignment != BondStereopermutator::Alignment::BetweenEclipsedAndStaggered
    ) {
      alignment = BondStereopermutator::Alignment::Eclipsed;
      if(bondType == BondType::Single) {
        alignment = BondStereopermutator::Alignment::Staggered;
      }
    }

    BondStereopermutator trialPermutator {
      *firstAtomStereopermutatorOption,
      *secondAtomStereopermutatorOption,
      bondIndex,
      alignment
    };

    if(trialPermutator.numAssignments() < 2) {
      return IgnoreReason::RotationIsIsotropic;
    }

    return trialPermutator;
  } catch(...) {
    return IgnoreReason::AtomStereopermutatorPreconditionsUnmet;
  }
}

DirectedConformerGenerator::Impl::Impl(
  Molecule molecule,
  const BondStereopermutator::Alignment alignment,
  const BondList& bondsToConsider
) : molecule_(std::move(molecule)), alignment_(alignment)
{
  // Allocate some space for the relevant bonds and new stereopermutators
  relevantBonds_.reserve(molecule_.graph().B() / 2);
  std::vector<BondStereopermutator> bondStereopermutators;
  bondStereopermutators.reserve(molecule_.graph().B() / 2);

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

  const auto processBond = [&](const BondIndex& bond) {
    auto importanceVariant = considerBond(bond, molecule_, alignment_);

    if(boost::apply_visitor(visitor, importanceVariant)) {
      relevantBonds_.push_back(bond);
    }
  };

  if(bondsToConsider.empty()) {
    Temple::forEach(molecule_.graph().bonds(), processBond);
  } else {
    Temple::forEach(bondsToConsider, processBond);
  }

  // Sort the relevant bonds and shrink
  Temple::sort(relevantBonds_);
  relevantBonds_.shrink_to_fit();

  /* In case there are no bonds to consider, then we're done. No other work
   * to be done.
   */
  if(relevantBonds_.empty()) {
    return;
  }

  // Add the BondStereopermutators to our underlying molecule's stereopermutator list
  for(auto& bondStereopermutator : bondStereopermutators) {
    molecule_.pImpl_->stereopermutators_.add(
      std::move(bondStereopermutator)
    );
  }

  // Collect the bounds on each stereopermutator's permutations for the trie
  auto bounds = Temple::map(
    relevantBonds_,
    [&](const BondIndex& bond) -> std::uint8_t {
      auto stereoOption = molecule_.stereopermutators().option(bond);
      if(!stereoOption) {
        throw std::logic_error(
          "No BondStereopermutator for relevant bond, but it was just inserted"
        );
      }

      return stereoOption->numAssignments();
    }
  );

  // Initialize the trie with the list of bounds
  decisionLists_.setBounds(std::move(bounds));
}

DirectedConformerGenerator::DecisionList
DirectedConformerGenerator::Impl::generateNewDecisionList(Random::Engine& engine) {
  if(relevantBonds_.empty()) {
    throw std::logic_error("List of relevant bonds is empty!");
  }

  Detail::BoundedNodeTrieChooseFunctor<std::uint8_t> chooseFunctor {engine};
  return decisionLists_.generateNewEntry(chooseFunctor);
}

Molecule DirectedConformerGenerator::Impl::conformationMolecule(const DecisionList& decisionList) const {
  StereopermutatorList permutators = molecule_.stereopermutators();

  const unsigned U = decisionList.size();
  if(U != decisionLists_.bounds().size() || U != relevantBonds_.size()) {
    throw std::invalid_argument("Passed decision list has wrong length");
  }

  Temple::forEach(
    Temple::Adaptors::zip(relevantBonds_, decisionList),
    [&](const BondIndex& bond, const std::uint8_t assignment) {
      permutators.option(bond)->assign(assignment);
    }
  );

  return Molecule(
    molecule_.graph(),
    std::move(permutators)
  );
}

outcome::result<Utils::PositionCollection>
DirectedConformerGenerator::Impl::checkGeneratedConformation(
  outcome::result<Utils::PositionCollection> conformerResult,
  const DecisionList& decisionList,
  const BondStereopermutator::FittingMode fitting
) const {
  if(conformerResult) {
    const auto generatedDecisionList = getDecisionList(
      conformerResult.value(),
      fitting
    );

    if(generatedDecisionList != decisionList) {
      return DgError::DecisionListMismatch;
    }
  }

  return conformerResult;
}

outcome::result<Utils::PositionCollection>
DirectedConformerGenerator::Impl::generateRandomConformation(
  const DecisionList& decisionList,
  const DistanceGeometry::Configuration& configuration,
  const BondStereopermutator::FittingMode fitting
) const {
  return checkGeneratedConformation(
    Scine::Molassembler::generateRandomConformation(
      conformationMolecule(decisionList),
      configuration
    ),
    decisionList,
    fitting
  );
}

outcome::result<Utils::PositionCollection>
DirectedConformerGenerator::Impl::generateConformation(
  const DecisionList& decisionList,
  const unsigned seed,
  const DistanceGeometry::Configuration& configuration,
  const BondStereopermutator::FittingMode fitting
) const {
  return checkGeneratedConformation(
    Scine::Molassembler::generateConformation(
      conformationMolecule(decisionList),
      seed,
      configuration
    ),
    decisionList,
    fitting
  );
}

DirectedConformerGenerator::DecisionList
DirectedConformerGenerator::Impl::getDecisionList(
  const Utils::AtomCollection& atomCollection,
  const BondStereopermutator::FittingMode fitting
) const {
  if(
    !Temple::all_of(
      molecule_.graph().atoms(),
      [&](const AtomIndex i) -> bool {
        return atomCollection.getElement(i) == molecule_.graph().elementType(i);
      }
    )
  ) {
    throw std::logic_error("Input AtomCollection elements do not match generator's underlying molecule. Misordered? Different molecule input?");
  }

  return getDecisionList(atomCollection.getPositions(), fitting);
}

DirectedConformerGenerator::DecisionList
DirectedConformerGenerator::Impl::getDecisionList(
  const Utils::PositionCollection& positions,
  const BondStereopermutator::FittingMode fitting
) const {
  const AngstromPositions angstromPositions {positions};

  if(
    !Temple::all_of(
      relevantBonds_,
      [&](const BondIndex& bondIndex) -> bool {
        const auto& permutators = molecule_.stereopermutators();
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

  /* Refit all atom stereopermutators and ensure stereopermutations are
   * identical, storing fitted shape maps for bond stereopermutator fitting later
   */
  std::unordered_map<AtomIndex, AtomStereopermutator::ShapeMap> shapeMaps;
  for(
    const AtomStereopermutator& stereopermutator :
    molecule_.stereopermutators().atomStereopermutators()
  ) {
    AtomStereopermutator refitted = stereopermutator;
    auto shapeMap = refitted.fit(molecule_.graph(), angstromPositions);
    if(refitted.getShape() != stereopermutator.getShape()) {
      const std::string error = (
        Shapes::name(refitted.getShape())
        + " was found instead of "
        + Shapes::name(stereopermutator.getShape())
        + " at atom "
        + std::to_string(stereopermutator.placement())
        + "! This indicates a precondition violation."
      );
      throw std::logic_error(error);
    }
    if(stereopermutator.assigned() && refitted.assigned() != stereopermutator.assigned()) {
      auto assignmentToString = [](const boost::optional<unsigned>& assignment) -> std::string {
        if(assignment) {
          return std::to_string(assignment.value());
        }

        return "U";
      };
      const std::string error = (
        "Assignment "
        + assignmentToString(refitted.assigned())
        + " was found instead of "
        + assignmentToString(stereopermutator.assigned())
        + " at atom "
        + std::to_string(stereopermutator.placement())
        + "! This indicates a precondition violation."
      );
      throw std::logic_error(error);
    }
    shapeMaps.emplace(stereopermutator.placement(), std::move(shapeMap.value()));
  }

  return Temple::map(
    relevantBonds_,
    [&](const BondIndex& bondIndex) -> std::uint8_t {
      const auto stereoOption = molecule_.stereopermutators().option(bondIndex);
      if(!stereoOption) {
        throw std::logic_error(
          "No BondStereopermutator for selected bond in molecule!"
        );
      }

      const auto fittingReferences = Temple::map(
        bondIndex,
        [&](const AtomIndex v) -> BondStereopermutator::FittingReferences {
          return {
            molecule_.stereopermutators().at(v),
            shapeMaps.at(v)
          };
        }
      );

      BondStereopermutator stereopermutator = stereoOption.value();
      stereopermutator.fit(angstromPositions, fittingReferences, fitting);
      return stereopermutator.assigned().value_or(unknownDecision);
    }
  );
}

void DirectedConformerGenerator::Impl::enumerate(
  std::function<void(const DecisionList&, Utils::PositionCollection)> callback,
  const unsigned seed,
  const EnumerationSettings& settings
) {
  clear();
  const unsigned size = idealEnsembleSize();

#pragma omp parallel for
  for(unsigned increment = 0; increment < size; ++increment) {
    Random::Engine localEngine(seed + increment);

    DecisionList decisionList;
#pragma omp critical(decisionSetAccess)
    {
      decisionList = generateNewDecisionList(localEngine);
    }

    for(unsigned i = 0; i < settings.dihedralRetries; ++i) {
      outcome::result<Utils::PositionCollection> conformer {DgError::DecisionListMismatch};

      try {
        conformer = generateConformation(
          decisionList,
          localEngine(),
          settings.configuration,
          settings.fitting
        );
      } catch(...) {}

      if(conformer) {
#pragma omp critical(guardCallback)
        {
          callback(decisionList, conformer.value());
        }
        break;
      }

      if(conformer.error() != DgError::DecisionListMismatch) {
        /* Only allow decision list failure retries for retries, break on
         * anything else
         */
        break;
      }
    }
  }
}

DirectedConformerGenerator::Relabeler DirectedConformerGenerator::Impl::relabeler() const {
  return Relabeler(relevantBonds_, molecule_);
}

std::vector<int>
DirectedConformerGenerator::Impl::binMidpointIntegers(
  const DecisionList& decisions
) const {
  return Temple::map(
    Temple::Adaptors::zip(relevantBonds_, decisions),
    [&](const BondIndex bond, const unsigned stereopermutation) -> int {
      const auto& permutator = molecule_.stereopermutators().at(bond);
      const auto& dominantDihedral = permutator.composite().allPermutations().at(stereopermutation).dihedrals.front();
      return static_cast<int>(std::round(Temple::Math::toDegrees(std::get<2>(dominantDihedral))));
    }
  );
}

std::vector<std::pair<int, int>>
DirectedConformerGenerator::Impl::binBounds(
  const DecisionList& decisions
) const {
  return Temple::map(
    Temple::Adaptors::zip(relevantBonds_, decisions),
    [&](const BondIndex bond, const unsigned stereopermutation) -> std::pair<int, int> {
      const auto& permutator = molecule_.stereopermutators().at(bond);
      const unsigned numStereopermutations = permutator.numStereopermutations();

      std::array<unsigned, 3> permutations {{
        (stereopermutation + (numStereopermutations - 1)) % numStereopermutations,
        stereopermutation,
        (stereopermutation + 1) % numStereopermutations
      }};

      const auto dominantAngles = Temple::map(permutations, [&](const unsigned permutation) -> double {
        const auto& dominantDihedral = permutator.composite().allPermutations().at(permutation).dihedrals.front();
        return std::get<2>(dominantDihedral);
      });

      const auto averages = Temple::map(
        Temple::Adaptors::sequentialPairs(dominantAngles),
        [](const double x, const double y) -> int {
          const double average = [&]() {
            if(x <= y) {
              return (x + y) / 2;
            }

            return (x + y) / 2 + M_PI;
          }();

          return static_cast<int>(std::round(Temple::Math::toDegrees(average)));
        }
      );

      // Increment the front to make exclusive bin bounds
      const unsigned front = [](int x) -> int {
        if(x == 180) {
          return -179;
        }

        return x + 1;
      }(averages.front());

      return std::make_pair(front, averages.back());
    }
  );
}

} // namespace Molassembler
} // namespace Scine
