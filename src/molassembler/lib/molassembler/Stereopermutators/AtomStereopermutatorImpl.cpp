/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/AtomStereopermutatorImpl.h"

#include "boost/range/join.hpp"
#include "chemical_symmetries/Properties.h"
#include "chemical_symmetries/DynamicProperties.h"
#include "CyclicPolygons.h"
#include <Eigen/Dense>
#include "stereopermutation/GenerateUniques.h"
#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"
#include "temple/Optionals.h"
#include "temple/Random.h"
#include "temple/TinySet.h"
#include "temple/constexpr/Math.h"
#include "temple/constexpr/Numeric.h"

#include "molassembler/Cycles.h"
#include "molassembler/Detail/BuildTypeSwitch.h"
#include "molassembler/Detail/DelibHelpers.h"
#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Log.h"
#include "molassembler/Modeling/CommonTrig.h"

namespace Scine {

namespace molassembler {

/* Static functions */

namespace detail {

Symmetry::Name pickTransition(
  const Symmetry::Name symmetryName,
  const unsigned T,
  boost::optional<unsigned> removedSymmetryPositionOptional
) {
  boost::optional<Symmetry::properties::SymmetryTransitionGroup> bestTransition;
  std::vector<Symmetry::Name> propositions;

  auto replaceOrAdd = [&](
    Symmetry::Name newName,
    const Symmetry::properties::SymmetryTransitionGroup& transitionGroup
  ) {
    if(bestTransition){
      if(transitionGroup.angularDistortion < bestTransition->angularDistortion) {
        bestTransition = transitionGroup;
        propositions = {newName};
      } else if(transitionGroup.angularDistortion == bestTransition->angularDistortion) {
        propositions.push_back(newName);
      }
    }

    if(!bestTransition) {
      bestTransition = transitionGroup;
      propositions = {newName};
    }
  };

  // Populate the list of symmetries with candidates
  for(const Symmetry::Name propositionalName : Symmetry::allNames) {
    if(Symmetry::size(propositionalName) != T) {
      continue;
    }

    if(auto transitionOptional = Symmetry::getMapping(symmetryName, propositionalName, removedSymmetryPositionOptional)) {
      if(transitionOptional->indexMappings.empty()) {
        continue;
      }

      if(
        Options::chiralStatePreservation == ChiralStatePreservation::EffortlessAndUnique
        && transitionOptional->indexMappings.size() == 1
        && transitionOptional->angularDistortion <= 0.2
      ) {
        replaceOrAdd(propositionalName, *transitionOptional);
      }

      if(
        Options::chiralStatePreservation == ChiralStatePreservation::Unique
        && transitionOptional->indexMappings.size() == 1
      ) {
        replaceOrAdd(propositionalName, *transitionOptional);
      }

      if(Options::chiralStatePreservation == ChiralStatePreservation::RandomFromMultipleBest) {
        replaceOrAdd(propositionalName, *transitionOptional);
      }
    }
  }

  /* In case we have no propositions, add all of fitting size. We need to
   * propose something, even if state will not be propagated due to the
   * preservation criteria.
   */
  if(propositions.empty()) {
    return Symmetry::properties::mostSymmetric(T);
  }

  return Symmetry::properties::mostSymmetric(std::move(propositions));
}

} // namespace detail

Symmetry::Name AtomStereopermutator::Impl::up(const Symmetry::Name symmetryName) {
  return detail::pickTransition(symmetryName, Symmetry::size(symmetryName) + 1, boost::none);
}

Symmetry::Name AtomStereopermutator::Impl::down(const Symmetry::Name symmetryName, const unsigned removedSymmetryPosition) {
  return detail::pickTransition(symmetryName, Symmetry::size(symmetryName) - 1, removedSymmetryPosition);
}

/* Constructors */
AtomStereopermutator::Impl::Impl(
  const OuterGraph& graph,
  // The symmetry of this Stereopermutator
  const Symmetry::Name symmetry,
  // The atom this Stereopermutator is centered on
  const AtomIndex centerAtom,
  // Ranking information of substituents
  RankingInformation ranking
) : _ranking {std::move(ranking)},
    _centerAtom {centerAtom},
    _symmetry {symmetry},
    _assignmentOption {boost::none}
{
  _cache = PermutationState {
    _ranking,
    _centerAtom,
    _symmetry,
    graph
  };
}

/* Modification */
void AtomStereopermutator::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment) {
    assert(assignment.value() < _cache.feasiblePermutations.size());
  }

  // Store current assignment
  _assignmentOption = std::move(assignment);

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndex -> unsigned).
   */
  if(_assignmentOption) {
    _cache.symmetryPositionMap = PermutationState::generateLigandToSymmetryPositionMap(
      _cache.permutations.stereopermutations.at(
        _cache.feasiblePermutations.at(
          _assignmentOption.value()
        )
      ),
      _cache.canonicalLigands
    );
  } else { // Wipe the map
    _cache.symmetryPositionMap.clear();
  }
}

void AtomStereopermutator::Impl::assignRandom() {
  assign(
    temple::random::pickDiscrete(
      // Map the feasible permutations onto their weights
      temple::map(
        _cache.feasiblePermutations,
        [&](const unsigned permutationIndex) -> unsigned {
          return _cache.permutations.weights.at(permutationIndex);
        }
      ),
      randomnessEngine()
    )
  );
}

void AtomStereopermutator::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // RankingInformation changes (lots of atom indices)
  _ranking.applyPermutation(permutation);

  // _centerAtom must change
  _centerAtom = permutation.at(_centerAtom);

  // Neither symmetry nor assignment change

  /* Although ranking and central atom are implicated in its creation,
   * PermutationState's state is completely independent of atom indices.
   */
}

boost::optional<AtomStereopermutator::PropagatedState> AtomStereopermutator::Impl::propagate(
  const OuterGraph& graph,
  RankingInformation newRanking,
  boost::optional<Symmetry::Name> symmetryOption
) {
  // If nothing changes, nothing changes.
  if(newRanking == _ranking) {
    return boost::none;
  }

#ifndef NDEBUG
  /* The following is a precondition. It is also rankPriority's postcondition,
   * but propagation may mess with it, so we need to be sure.
   */
  for(const auto& ligand : boost::range::join(_ranking.ligands, newRanking.ligands)) {
    assert(
      std::is_sorted(
        std::begin(ligand),
        std::end(ligand)
      )
    );
  }
#endif

  /* There are five possible situations that can occur here:
   * - A substituent is removed, and that substituent was sole constituent of a
   *   ligand, leading to a symmetry size decrease.
   * - A substituent is removed, but that substituent was part of a ligand with
   *   multiple constituents. The symmetry stays the same.
   * - No substituent is removed or added, but there is a ranking change. The
   *   symmetry stays the same.
   * - A substituent is added, but to an existing ligand. The symmetry stays
   *   the same.
   * - A substituent is added that solely constitutes a new ligand. The
   *   symmetry size increases.
   */
  enum class PropagationSituation {
    Unknown,
    LigandRemoval,
    SubstituentRemoval,
    RankingChange,
    SubstituentAddition,
    LigandAddition
  };

  /* Fully determine in which situation we are and all conditional information
   * needed to proceed with each.
   */
  auto situation = PropagationSituation::Unknown;
  int ligandCountChange = static_cast<int>(newRanking.ligands.size() - _ranking.ligands.size());
  boost::optional<unsigned> alteredLigandIndex;
  boost::optional<AtomIndex> alteredSubstituentIndex;

  // Complain if more than one ligand is changed either way
  if(std::abs(ligandCountChange) > 1) {
    throw std::logic_error("propagateGraphChange should only ever handle a single ligand addition or removal at once");
  }

  auto countSubstituents = [](const RankingInformation::LigandsType& ligands) -> unsigned {
    return temple::accumulate(
      ligands,
      0u,
      [](const unsigned carry, const auto& ligand) -> unsigned {
        return carry + ligand.size();
      }
    );
  };

  int substituentCountChange = static_cast<int>(countSubstituents(newRanking.ligands)) - static_cast<int>(countSubstituents(_ranking.ligands));
  // Complain if multiple ligands are changed either way
  if(std::abs(substituentCountChange) > 1) {
    throw std::logic_error("propagateGraphChange should only ever handle a single subtituent addition or removal at once");
  }

  // Collect all substituents of a ranking and sort them
  auto collectSubstituents = [](const RankingInformation::LigandsType& ligands) -> std::vector<AtomIndex> {
    std::vector<AtomIndex> substituents;
    substituents.reserve(2 * ligands.size());

    for(const auto& ligand : ligands) {
      for(const AtomIndex& index : ligand) {
        substituents.push_back(index);
      }
    }

    substituents.shrink_to_fit();

    temple::inplace::sort(substituents);

    return substituents;
  };

  auto determineChangedSubstituentAndLigand = [&](
    const RankingInformation& rankingWithMoreSubstituents,
    const RankingInformation& rankingWithFewerSubstituents
  ) {
    // Figure out which substituent has changed, if any
    auto largerSubstituentList = collectSubstituents(rankingWithMoreSubstituents.ligands);
    auto smallerSubstituentList = collectSubstituents(rankingWithFewerSubstituents.ligands);

    std::vector<AtomIndex> changedSubstituents;

    std::set_difference(
      std::begin(largerSubstituentList),
      std::end(largerSubstituentList),
      std::begin(smallerSubstituentList),
      std::end(smallerSubstituentList),
      std::back_inserter(changedSubstituents)
    );

    /* Just because the count difference is +1 does not mean that only one has
     * changed, so we diligently recheck:
     */
    if(changedSubstituents.size() > 1) {
      throw std::logic_error("propagateGraphChange should only ever handle a single substituent addition or removal at once");
    }

    assert(!changedSubstituents.empty());

    alteredSubstituentIndex = changedSubstituents.front();
    alteredLigandIndex = rankingWithMoreSubstituents.getLigandIndexOf(changedSubstituents.front());
  };

  if(substituentCountChange == +1) {
    determineChangedSubstituentAndLigand(newRanking, _ranking);

    assert(!newRanking.ligands.at(*alteredLigandIndex).empty());
    if(newRanking.ligands.at(*alteredLigandIndex).size() > 1) {
      situation = PropagationSituation::SubstituentAddition;
    } else {
      situation = PropagationSituation::LigandAddition;
    }
  } else if(substituentCountChange == -1) {
    determineChangedSubstituentAndLigand(_ranking, newRanking);
    assert(!_ranking.ligands.at(*alteredLigandIndex).empty());
    if(_ranking.ligands.at(*alteredLigandIndex).size() > 1) {
      situation = PropagationSituation::SubstituentRemoval;
    } else {
      situation = PropagationSituation::LigandRemoval;
    }
  } else {
    situation = PropagationSituation::RankingChange;
  }

  /* Decide the new symmetry */
  Symmetry::Name newSymmetry = symmetryOption.value_or_eval(
    [&]() {
      if(ligandCountChange == +1) {
        return up(_symmetry);
      }

      if(ligandCountChange == 0) {
        return _symmetry;
      }

      assert(ligandCountChange == -1);
      assert(alteredLigandIndex);

      /* We can only figure out a mapping if the stereocenter is assigned
       * since otherwise _cache.symmetryPositionMap is empty and the symmetry
       * position being removed is a necessary argument to down.
       */
      if(_assignmentOption) {
        return down(
          _symmetry,
          _cache.symmetryPositionMap.at(alteredLigandIndex.value())
        );
      }

      // Just return the most symmetric symmetry of the target size
      return Symmetry::properties::mostSymmetric(
        newRanking.ligands.size()
      );
    }
  );

  /* Generate new assignments */
  PermutationState newPermutationState {
    newRanking,
    _centerAtom,
    newSymmetry,
    graph
  };

  boost::optional<unsigned> newStereopermutation = boost::none;

  /* Now we will attempt to propagate our assignment to the new set of
   * assignments. This is only necessary in the case that the stereopermutator is currently
   * assigned.
   *
   * In the case of an equal number of ligands, it is only possible to
   * propagate chiral state if the new number of assignments is smaller or
   * equal to the amount we have currently. This is because say we have AABCDE
   * in octahedral, we carry chiral state, and a ranking change leads to
   * ABCDEF.  There will be multiple ways to split A to F!
   */
  if(_assignmentOption && numStereopermutations() > 1) {
    // A flat map of symmetry position to new ligand index in the new symmetry
    std::vector<unsigned> ligandsAtNewSymmetryPositions;

    /* Ligand-level changes */
    if(situation == PropagationSituation::LigandAddition) {
      /* Try to get a mapping to the new symmetry. If that returns a Some, try
       * to get a mapping by preservationOption policy. If any of these steps
       * returns boost::none, the whole expression is boost::none.
       */
      auto suitableMappingOption = Symmetry::getMapping(
        _symmetry,
        newSymmetry,
        boost::none
      ) | temple::callIfSome(
        PermutationState::getIndexMapping,
        temple::ANS,
        Options::chiralStatePreservation
      );

      if(suitableMappingOption) {
        /* So now we must transfer the current assignment into the new symmetry
         * and search for it in the set of uniques.
         */
        const auto& symmetryMapping = suitableMappingOption.value();

        // Apply the mapping to get old symmetry positions at their places in the new symmetry
        auto oldSymmetryPositionsInNewSymmetry = Symmetry::properties::applyIndexMapping(
          newSymmetry,
          symmetryMapping
        );

        // Invert symmetryPositionMap to get: ligand = map.at(symmetryPosition)
        std::vector<unsigned> oldSymmetryPositionToLigandMap(_cache.symmetryPositionMap.size());
        for(unsigned i = 0 ; i < _cache.symmetryPositionMap.size(); ++i) {
          oldSymmetryPositionToLigandMap.at(
            _cache.symmetryPositionMap.at(i)
          ) = i;
        }
        /* We assume the new ligand is also merely added to the end! This may
         * not always be true
         */
        oldSymmetryPositionToLigandMap.push_back(_cache.symmetryPositionMap.size());

        // Replace old symmetry positions by their ligand indices
        ligandsAtNewSymmetryPositions = temple::map(
          oldSymmetryPositionsInNewSymmetry,
          [&oldSymmetryPositionToLigandMap](const unsigned oldSymmetryPosition) -> unsigned {
            return oldSymmetryPositionToLigandMap.at(oldSymmetryPosition);
          }
        );
      }
      /* If no mapping can be found that fits to the preservationOption,
       * newStereopermutation remains boost::none, and this stereopermutator loses
       * any chiral information it may have had.
       */
    }

    if(situation == PropagationSituation::LigandRemoval) {
      /* Try to get a mapping to the new symmetry. If that returns a Some, try
       * to get a mapping by preservationOption policy. If any of these steps
       * returns boost::none, the whole expression is boost::none.
       */
      auto suitableMappingOptional = Symmetry::getMapping(
        _symmetry,
        newSymmetry,
        /* Last parameter is the deleted symmetry position, which is the
         * symmetry position at which the ligand being removed is currently at
         */
        _cache.symmetryPositionMap.at(*alteredLigandIndex)
      ) | temple::callIfSome(
        PermutationState::getIndexMapping,
        temple::ANS,
        Options::chiralStatePreservation
      );

      /* It is weird that this block of code works although it is structurally
       * so different from ligand addition propagation. Is it coincidence?
       */
      if(suitableMappingOptional) {
        /* symmetryMapping {0, 4, 3, 2} means
         * - at symmetry position 0 in the new symmetry, we have what was
         *   at symmetry position 0 in the old symmetry
         * - at symmetry position 1 in the new symmetry, we have what was
         *   at symmetry position 4 in the old symmetry
         *
         */
        const auto& symmetryMapping = suitableMappingOptional.value();

        // Invert symmetryPositionMap to get: ligand = map.at(symmetryPosition)
        std::vector<unsigned> oldSymmetryPositionToLigandMap(_cache.symmetryPositionMap.size());
        for(unsigned i = 0 ; i < _cache.symmetryPositionMap.size(); ++i) {
          oldSymmetryPositionToLigandMap.at(
            _cache.symmetryPositionMap.at(i)
          ) = i;
        }

        // Transfer indices from current symmetry to new symmetry
        ligandsAtNewSymmetryPositions.resize(Symmetry::size(newSymmetry));
        for(unsigned i = 0; i < Symmetry::size(newSymmetry); ++i) {
          // i is a symmetry position
          ligandsAtNewSymmetryPositions.at(i) = oldSymmetryPositionToLigandMap.at(
            symmetryMapping.at(i)
          );
        }
      }
    }

    /* Substituent-level changes */
    if(situation == PropagationSituation::SubstituentAddition) {
      /* Sort ligands' constituting atom indices in both rankings so we can use
       * lexicographical comparison to create a mapping. Add the
       * newSubstituentIndex to the old ranking's appropriate ligand so the
       * mapping is 1:1.
       */
      temple::TinySet<AtomIndex>::checked_insert(
        _ranking.ligands.at(*alteredLigandIndex),
        *alteredSubstituentIndex
      );

      // Generate a flat map of old ligand indices to new ligand indices
      auto ligandMapping = temple::map(
        _ranking.ligands,
        [&newRanking](const auto& ligand) -> unsigned {
          auto findIter = std::find(
            newRanking.ligands.begin(),
            newRanking.ligands.end(),
            ligand
          );

          assert(findIter != newRanking.ligands.end());

          return findIter - newRanking.ligands.begin();
        }
      );

      // Write the found mapping into ligandsAtNewSymmetryPositions
      ligandsAtNewSymmetryPositions.resize(Symmetry::size(newSymmetry));
      for(unsigned i = 0; i < ligandMapping.size(); ++i) {
        ligandsAtNewSymmetryPositions.at(i) = ligandMapping.at(
          _cache.symmetryPositionMap.at(i)
        );
      }

      /* Revert the change to _ranking from earlier for easy mapping so that
       * the returned prior state is unchanged.
       */
      temple::inplace::remove(
        _ranking.ligands.at(*alteredLigandIndex),
        *alteredSubstituentIndex
      );
    }

    if(situation == PropagationSituation::SubstituentRemoval) {
      /* Sort ligands in the old ranking and new so we can use lexicographical
       * comparison to figure out a mapping
       */
      temple::inplace::remove(
        _ranking.ligands.at(*alteredLigandIndex),
        *alteredSubstituentIndex
      );

      // Calculate the mapping from old ligands to new ones
      auto ligandMapping = temple::map(
        _ranking.ligands,
        [&newRanking](const auto& ligand) -> unsigned {
          auto findIter = std::find(
            newRanking.ligands.begin(),
            newRanking.ligands.end(),
            ligand
          );

          assert(findIter != newRanking.ligands.end());

          return findIter - newRanking.ligands.begin();
        }
      );

      ligandsAtNewSymmetryPositions.resize(Symmetry::size(newSymmetry));
      // Transfer ligands to new mapping
      for(unsigned i = 0; i < ligandMapping.size(); ++i) {
        ligandsAtNewSymmetryPositions.at(i) = ligandMapping.at(
          _cache.symmetryPositionMap.at(i)
        );
      }

      /* Revert the change to _ranking from earlier for easy mapping so that
       * the returned prior state is unchanged.
       */
      temple::TinySet<AtomIndex>::checked_insert(
        _ranking.ligands.at(*alteredLigandIndex),
        *alteredSubstituentIndex
      );
    }

    /* Ranking-level change */
    if(
      situation == PropagationSituation::RankingChange
      && (
        newPermutationState.permutations.stereopermutations.size()
        <= _cache.permutations.stereopermutations.size()
      )
    ) {
      const auto& currentStereopermutation = _cache.permutations.stereopermutations.at(
        _cache.feasiblePermutations.at(
          _assignmentOption.value()
        )
      );

      // Replace the characters by their corresponding indices from the old ranking
      ligandsAtNewSymmetryPositions = PermutationState::generateSymmetryPositionToLigandMap(
        currentStereopermutation,
        _cache.canonicalLigands
      );
    }

    if(!ligandsAtNewSymmetryPositions.empty()) {
      // Replace the ligand indices by their new ranking characters
      auto newStereopermutationCharacters = PermutationState::makeStereopermutationCharacters(
        newPermutationState.canonicalLigands,
        newPermutationState.symbolicCharacters,
        ligandsAtNewSymmetryPositions
      );

      // Create a new assignment with those characters
      auto trialStereopermutation = stereopermutation::Stereopermutation(
        newStereopermutationCharacters,
        newPermutationState.selfReferentialLinks
      );

      // Generate all rotations of this trial assignment
      auto allTrialRotations = trialStereopermutation.generateAllRotations(newSymmetry);

      // Find out which of the new assignments has a rotational equivalent
      for(unsigned i = 0; i < newPermutationState.permutations.stereopermutations.size(); ++i) {
        if(allTrialRotations.count(newPermutationState.permutations.stereopermutations.at(i)) > 0) {
          newStereopermutation = i;
          break;
        }
      }
    }
  }

  auto oldStateTuple = std::make_tuple(
    std::move(_ranking),
    std::move(_cache),
    std::move(_assignmentOption)
  );

  // Overwrite the class state
  _ranking = std::move(newRanking);
  _symmetry = newSymmetry;
  _cache = std::move(newPermutationState);
  assign(newStereopermutation);

  return {std::move(oldStateTuple)};
}

void AtomStereopermutator::Impl::propagateVertexRemoval(const AtomIndex removedIndex) {
  /* This function replaces any occurrences of the atom index that is being
   * removed in the global state with a placeholder of the same type and updates
   * any invalidated atom indices.
   */

  /* If the central atom is being removed, just drop this stereopermutator
   * beforehand in caller. This would just be unnecessary work.
   */
  assert(_centerAtom != removedIndex);

  // Define some helper functions
  auto updateIndexInplace = [&removedIndex](AtomIndex& index) -> void {
    if(index > removedIndex) {
      --index;
    } else if(index == removedIndex) {
      index = InnerGraph::removalPlaceholder;
    }
  };

  auto updateIndex = [&removedIndex](const AtomIndex index) -> AtomIndex {
    if(index > removedIndex) {
      return index - 1;
    }

    if(index == removedIndex) {
      return InnerGraph::removalPlaceholder;
    }

    return index;
  };

  /* Update indices in RankingInformation */
  for(auto& equalPrioritySet : _ranking.sortedSubstituents) {
    for(auto& index : equalPrioritySet) {
      updateIndexInplace(index);
    }
  }

  for(auto& ligandIndicesList : _ranking.ligands) {
    for(auto& atomIndex : ligandIndicesList) {
      updateIndexInplace(atomIndex);
    }
  }

  for(auto& link : _ranking.links) {
    link.cycleSequence = temple::map(
      link.cycleSequence,
      updateIndex
    );
  }
}

const PermutationState& AtomStereopermutator::Impl::getPermutationState() const {
  return _cache;
}

const RankingInformation& AtomStereopermutator::Impl::getRanking() const {
  return _ranking;
}

Symmetry::Name AtomStereopermutator::Impl::getSymmetry() const {
  return _symmetry;
}

std::vector<unsigned> AtomStereopermutator::Impl::getSymmetryPositionMap() const {
  if(_assignmentOption == boost::none) {
    throw std::logic_error(
      "The AtomStereopermutator is unassigned, ligands are not assigned to "
      "symmetry positions"
    );
  }

  return _cache.symmetryPositionMap;
}

void AtomStereopermutator::Impl::fit(
  const OuterGraph& graph,
  const AngstromWrapper& angstromWrapper
) {
  const unsigned S = Symmetry::size(_symmetry);

  // For all atoms making up a ligand, decide on the spatial average position
  const std::vector<Eigen::Vector3d> ligandPositions = temple::map(
    _ranking.ligands,
    [&angstromWrapper](const std::vector<AtomIndex>& ligandAtoms) -> Eigen::Vector3d {
      return DelibHelpers::averagePosition(angstromWrapper.positions, ligandAtoms);
    }
  );

  std::vector<Symmetry::Name> excludeSymmetries;

  /* Special case for trigonal bipyramidal and square pyramidal: tau-calculation
   */
  if(
    Options::tauCriterion == TauCriterion::Enable
    && (S == 4 || S == 5)
  ) {
    // Find the two largest angles
    std::vector<double> angles;
    for(unsigned i = 0; i < S; ++i) {
      for(unsigned j = i + 1; j < S; ++j) {
        angles.push_back(
          DelibHelpers::angle(
            ligandPositions.at(i),
            angstromWrapper.positions.row(_centerAtom),
            ligandPositions.at(j)
          )
        );
      }
    }

    std::sort(std::begin(angles), std::end(angles));

    // Largest angle is beta
    const double beta = angles.back();
    // Second-largest angle is alpha
    const double alpha = angles.at(angles.size() - 2);
    if(S == 4) {
      constexpr double theta = temple::Math::toRadians(109.5);
      const double tau_prime = (
        (beta - alpha) / (2 * M_PI - theta)
        + (M_PI - beta) / (M_PI - theta)
      );

      /* Significance
       * - τ₄' = 0 -> Symmetry is square planar
       * - τ₄' = 0.24 -> Symmetry is seesaw
       * - τ₄' = 1 -> Symmetry is tetrahedral
       */
      if(tau_prime < 0.12) {
        // Symmetry is square planar
        excludeSymmetries.push_back(Symmetry::Name::Seesaw);
        excludeSymmetries.push_back(Symmetry::Name::Tetrahedral);
      } else if(0.12 <= tau_prime && tau_prime < 0.62) {
        excludeSymmetries.push_back(Symmetry::Name::SquarePlanar);
        // Symmetry is seesaw
        excludeSymmetries.push_back(Symmetry::Name::Tetrahedral);
      } else if(0.62 <= tau_prime) {
        excludeSymmetries.push_back(Symmetry::Name::SquarePlanar);
        excludeSymmetries.push_back(Symmetry::Name::Seesaw);
        // Symmetry is tetrahedral
      }
    } else if(S == 5) {
      // Calculate tau as largest angle minus second-largest angle divided by 60°
      const double tau = (angles.back() - angles.at(angles.size() - 2)) / temple::Math::toRadians(60.0);

      /* Significance:
       * - τ₅ = 0 -> Symmetry is square pyramidal
       * - τ₅ = 1 -> Symmetry is trigonal bipyramidal
       */

      if(tau < 0.5) {
        excludeSymmetries.push_back(Symmetry::Name::TrigonalBiPyramidal);
      } else if(tau > 0.5) {
        excludeSymmetries.push_back(Symmetry::Name::SquarePyramidal);
      }
    }
  }

  // Save stereopermutator state to return to if no fit is viable
  const Symmetry::Name priorSymmetry = _symmetry;
  const boost::optional<unsigned> priorStereopermutation  = _assignmentOption;

  const Symmetry::Name initialSymmetry {Symmetry::Name::Linear};
  const unsigned initialStereopermutation = 0;
  const double initialPenalty = 100;

  Symmetry::Name bestSymmetry = initialSymmetry;
  unsigned bestStereopermutation = initialStereopermutation;
  double bestPenalty = initialPenalty;
  unsigned bestStereopermutationMultiplicity = 1;

  auto excludesContains = temple::makeContainsPredicate(excludeSymmetries);

  // Cycle through all symmetries
  for(const auto& symmetryName : Symmetry::allNames) {
    // Skip any Symmetries of different size
    if(Symmetry::size(symmetryName) != S || excludesContains(symmetryName)) {
      continue;
    }

    // Change the symmetry of the AtomStereopermutator
    setSymmetry(symmetryName, graph);

    for(
      unsigned assignment = 0;
      assignment < numAssignments();
      ++assignment
    ) {
      // Assign the stereopermutator
      assign(assignment);

      const double angleDeviations = temple::sum(
        temple::adaptors::transform(
          temple::adaptors::allPairs(
            temple::adaptors::range(Symmetry::size(_symmetry))
          ),
          [&](const unsigned ligandI, const unsigned ligandJ) -> double {
            return std::fabs(
              DelibHelpers::angle(
                ligandPositions.at(ligandI),
                angstromWrapper.positions.row(_centerAtom),
                ligandPositions.at(ligandJ)
              ) - angle(ligandI, ligandJ)
            );
          }
        )
      );

      // We can stop immediately if this is worse
      if(angleDeviations > bestPenalty) {
        continue;
      }

      /*! @todo should this be kept at all? Just a follow-up error from the angle
       * What value does it bring?
       */
      const double oneThreeDistanceDeviations = temple::sum(
        temple::adaptors::transform(
          temple::adaptors::allPairs(
            temple::adaptors::range(Symmetry::size(_symmetry))
          ),
          [&](const unsigned ligandI, const unsigned ligandJ) -> double {
            return std::fabs(
              // ligandI - ligandJ 1-3 distance from positions
              DelibHelpers::distance(
                ligandPositions.at(ligandI),
                ligandPositions.at(ligandJ)
              )
              // idealized 1-3 distance from
              - CommonTrig::lawOfCosines(
                // i-j 1-2 distance from positions
                DelibHelpers::distance(
                  ligandPositions.at(ligandI),
                  angstromWrapper.positions.row(_centerAtom)
                ),
                // j-k 1-2 distance from positions
                DelibHelpers::distance(
                  angstromWrapper.positions.row(_centerAtom),
                  ligandPositions.at(ligandJ)
                ),
                // idealized Stereopermutator angle
                angle(ligandI, ligandJ)
              )
            );
          }
        )
      );

      // Another early continue
      if(angleDeviations + oneThreeDistanceDeviations > bestPenalty) {
        continue;
      }

      const double chiralityDeviations = temple::sum(
        temple::adaptors::transform(
          minimalChiralityConstraints(),
          [&](const auto& minimalPrototype) -> double {
            auto fetchPosition = [&](const boost::optional<unsigned>& ligandIndexOptional) -> Eigen::Vector3d {
              if(ligandIndexOptional) {
                return ligandPositions.at(ligandIndexOptional.value());
              }

              return angstromWrapper.positions.row(_centerAtom);
            };

            double volume = DelibHelpers::adjustedSignedVolume(
              fetchPosition(minimalPrototype[0]),
              fetchPosition(minimalPrototype[1]),
              fetchPosition(minimalPrototype[2]),
              fetchPosition(minimalPrototype[3])
            );

            // minimalChiralityConstraints() supplies only Positive targets
            if(volume < 0) {
              return 1;
            }

            return 0;
          }
        )
      );

      double fitPenalty = angleDeviations
        + oneThreeDistanceDeviations
        + chiralityDeviations;


#ifndef NDEBUG
      Log::log(Log::Particulars::AtomStereopermutatorFit)
        << Symmetry::nameIndex(symmetryName)
        << ", " << assignment
        << ", " << std::setprecision(4) << std::fixed
        << angleDeviations << ", "
        << oneThreeDistanceDeviations << ", "
        << chiralityDeviations
        << std::endl;
#endif

      if(fitPenalty < bestPenalty) {
        bestSymmetry = symmetryName;
        bestStereopermutation = assignment;
        bestPenalty = fitPenalty;
        bestStereopermutationMultiplicity = 1;
      } else if(fitPenalty == bestPenalty) {
        // Assume that IF we have multiplicity, it's from the same symmetry
        assert(bestSymmetry == symmetryName);
        bestStereopermutationMultiplicity += 1;
      }
    }
  }

  /* In case NO assignments could be tested, return to the prior state.
   * This guards against situations in which predicates in
   * uniques could lead no assignments to be returned, such as
   * in e.g. square-planar AAAB with {0, 3}, {1, 3}, {2, 3} with removal of
   * trans-spanning groups. In that situation, all possible assignments are
   * trans-spanning and uniques is an empty vector.
   *
   * At the moment, this predicate is disabled, so no such issues should arise.
   * Just being safe.
   */
  if(
    bestSymmetry == initialSymmetry
    && bestStereopermutation == initialStereopermutation
    && bestPenalty == initialPenalty
  ) {
    // Return to prior
    setSymmetry(priorSymmetry, graph);
    assign(priorStereopermutation);
  } else {
    // Set to best fit
    setSymmetry(bestSymmetry, graph);

    /* How to handle multiplicity?
     * Current policy: If there is multiplicity, do not assign
     */
    if(bestStereopermutationMultiplicity > 1) {
      assign(boost::none);
    } else {
      assign(bestStereopermutation);
    }
  }
}

/* Information */
double AtomStereopermutator::Impl::angle(
  const unsigned i,
  const unsigned j
) const {
  assert(i != j);
  assert(!_cache.symmetryPositionMap.empty());

  return Symmetry::angleFunction(_symmetry)(
    _cache.symmetryPositionMap.at(i),
    _cache.symmetryPositionMap.at(j)
  );
}

boost::optional<unsigned> AtomStereopermutator::Impl::assigned() const {
  return _assignmentOption;
}

AtomIndex AtomStereopermutator::Impl::centralIndex() const {
  return _centerAtom;
}

boost::optional<unsigned> AtomStereopermutator::Impl::indexOfPermutation() const {
  if(_assignmentOption) {
    return _cache.feasiblePermutations.at(_assignmentOption.value());
  }

  return boost::none;
}

std::vector<
  std::array<boost::optional<unsigned>, 4>
> AtomStereopermutator::Impl::minimalChiralityConstraints() const {
  std::vector<
    std::array<boost::optional<unsigned>, 4>
  > precursors;

  // Only collect constraints if it's actually assigned
  if(_assignmentOption && numStereopermutations() > 1) {

    /* Invert _neighborSymmetryPositionMap, we need a mapping of
     *  (position in symmetry) -> atom index
     */
    auto symmetryPositionToLigandIndexMap = PermutationState::generateSymmetryPositionToLigandMap(
      _cache.permutations.stereopermutations.at(
        _cache.feasiblePermutations.at(
          _assignmentOption.value()
        )
      ),
      _cache.canonicalLigands
    );

    // Get list of tetrahedra from symmetry
    const auto& tetrahedraList = Symmetry::tetrahedra(_symmetry);

    precursors.reserve(tetrahedraList.size());
    for(const auto& tetrahedron : tetrahedraList) {
      /* Replace indices (represent positions within the symmetry) with the
       * ligand index at that position from the inverted map
       */

      // Make a minimal sequence from it
      precursors.push_back(
        temple::map(
          tetrahedron,
          [&](const boost::optional<unsigned>& indexOptional) -> boost::optional<unsigned> {
            if(indexOptional) {
              return symmetryPositionToLigandIndexMap.at(
                indexOptional.value()
              );
            }

            return boost::none;
          }
        )
      );
    }
  }

  return precursors;
}

std::string AtomStereopermutator::Impl::info() const {
  std::string returnString = "A on "s
    + std::to_string(_centerAtom) + " ("s + Symmetry::name(_symmetry) +", "s;

  const auto& characters = _cache.symbolicCharacters;
  std::copy(
    characters.begin(),
    characters.end(),
    std::back_inserter(returnString)
  );

  for(const auto& link : _cache.selfReferentialLinks) {
    returnString += ", "s + characters.at(link.first) + "-"s + characters.at(link.second);
  }

  returnString += "): "s;

  if(_assignmentOption) {
    returnString += std::to_string(_assignmentOption.value());
  } else {
    returnString += "u";
  }

  const unsigned A = numAssignments();
  returnString += "/"s + std::to_string(A);

  const unsigned P = numStereopermutations();
  if(P != A) {
    returnString += " ("s + std::to_string(P) + ")"s;
  }

  return returnString;
}

std::string AtomStereopermutator::Impl::rankInfo() const {
  /* rankInfo is specifically geared towards RankingTree's consumption,
   * and MUST use indices of permutation
   */
  return (
    "A-"s + std::to_string(static_cast<unsigned>(_symmetry))
    + "-"s + std::to_string(numStereopermutations())
    + "-"s + (
      indexOfPermutation()
      ? std::to_string(indexOfPermutation().value())
      : "u"s
    )
  );
}

unsigned AtomStereopermutator::Impl::numAssignments() const {
  return _cache.feasiblePermutations.size();
}

unsigned AtomStereopermutator::Impl::numStereopermutations() const {
  return _cache.permutations.stereopermutations.size();
}

void AtomStereopermutator::Impl::setSymmetry(
  const Symmetry::Name symmetryName,
  const OuterGraph& graph
) {
  if(_symmetry == symmetryName) {
    return;
  }

  _symmetry = symmetryName;

  _cache = PermutationState {
    _ranking,
    _centerAtom,
    _symmetry,
    graph
  };

  // Dis-assign the stereopermutator
  assign(boost::none);
}

bool AtomStereopermutator::Impl::operator == (const AtomStereopermutator::Impl& other) const {
  return (
    _symmetry == other._symmetry
    && _centerAtom == other._centerAtom
    && numStereopermutations() == other.numStereopermutations()
    && _assignmentOption == other._assignmentOption
  );
}

bool AtomStereopermutator::Impl::operator < (const AtomStereopermutator::Impl& other) const {
  unsigned thisPermutations = numStereopermutations(),
           otherPermutations = other.numStereopermutations();
  /* Sequentially compare individual components, comparing assignments last
   * if everything else matches
   */
  return (
    std::tie( _centerAtom, _symmetry, thisPermutations, _assignmentOption)
    < std::tie(other._centerAtom, other._symmetry, otherPermutations, other._assignmentOption)
  );
}

} // namespace molassembler

} // namespace Scine
