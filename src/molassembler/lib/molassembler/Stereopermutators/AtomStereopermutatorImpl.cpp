/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/AtomStereopermutatorImpl.h"

#include "boost/range/join.hpp"
#include "chemical_symmetries/DynamicProperties.h"
#include "chemical_symmetries/Properties.h"
#include "chemical_symmetries/TauCriteria.h"
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
#include "molassembler/Detail/Cartesian.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Log.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/Stereopermutators/ShapeVertexMaps.h"

namespace Scine {

namespace molassembler {

/* Static functions */

namespace detail {

Symmetry::Shape pickTransition(
  const Symmetry::Shape shape,
  const unsigned T,
  boost::optional<unsigned> removedSymmetryPositionOptional
) {
  boost::optional<Symmetry::properties::SymmetryTransitionGroup> bestTransition;
  std::vector<Symmetry::Shape> propositions;

  auto replaceOrAdd = [&](
    Symmetry::Shape newName,
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

  // Populate the list of shape with candidates
  for(const Symmetry::Shape propositionalShape : Symmetry::allShapes) {
    if(Symmetry::size(propositionalShape) != T) {
      continue;
    }

    if(auto transitionOptional = Symmetry::getMapping(shape, propositionalShape, removedSymmetryPositionOptional)) {
      if(transitionOptional->indexMappings.empty()) {
        continue;
      }

      if(
        Options::chiralStatePreservation == ChiralStatePreservation::EffortlessAndUnique
        && transitionOptional->indexMappings.size() == 1
        && transitionOptional->angularDistortion <= 0.2
      ) {
        replaceOrAdd(propositionalShape, *transitionOptional);
      }

      if(
        Options::chiralStatePreservation == ChiralStatePreservation::Unique
        && transitionOptional->indexMappings.size() == 1
      ) {
        replaceOrAdd(propositionalShape, *transitionOptional);
      }

      if(Options::chiralStatePreservation == ChiralStatePreservation::RandomFromMultipleBest) {
        replaceOrAdd(propositionalShape, *transitionOptional);
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

Symmetry::Shape AtomStereopermutator::Impl::up(const Symmetry::Shape shape) {
  return detail::pickTransition(shape, Symmetry::size(shape) + 1, boost::none);
}

Symmetry::Shape AtomStereopermutator::Impl::down(const Symmetry::Shape shape, const unsigned removedShapePosition) {
  return detail::pickTransition(shape, Symmetry::size(shape) - 1, removedShapePosition);
}

boost::optional<std::vector<unsigned>> AtomStereopermutator::Impl::getIndexMapping(
  const Symmetry::properties::SymmetryTransitionGroup& mappingsGroup,
  const ChiralStatePreservation& preservationOption
) {
  if(mappingsGroup.indexMappings.empty()) {
    return boost::none;
  }

  if(
    preservationOption == ChiralStatePreservation::EffortlessAndUnique
    && mappingsGroup.indexMappings.size() == 1
    && mappingsGroup.angularDistortion <= 0.2
  ) {
    return mappingsGroup.indexMappings.front();
  }

  if(
    preservationOption == ChiralStatePreservation::Unique
    && mappingsGroup.indexMappings.size() == 1
  ) {
    return mappingsGroup.indexMappings.front();
  }

  if(preservationOption == ChiralStatePreservation::RandomFromMultipleBest) {
    return mappingsGroup.indexMappings.at(
      temple::random::getSingle<unsigned>(
        0,
        mappingsGroup.indexMappings.size() - 1,
        randomnessEngine()
      )
    );
  }

  return boost::none;
}

/* Constructors */
AtomStereopermutator::Impl::Impl(
  const OuterGraph& graph,
  // The symmetry of this Stereopermutator
  const Symmetry::Shape shape,
  // The atom this Stereopermutator is centered on
  const AtomIndex centerAtom,
  // Ranking information of substituents
  RankingInformation ranking
) : _centerAtom {centerAtom},
    _shape {shape},
    _ranking {std::move(ranking)},
    _abstract {_ranking, _shape},
    _feasible {_abstract, _shape, _centerAtom, _ranking, graph},
    _assignmentOption {boost::none},
    _shapePositionMap {}
{}

/* Modification */
void AtomStereopermutator::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment) {
    assert(assignment.value() < _feasible.indices.size());
  }

  // Store new assignment
  _assignmentOption = std::move(assignment);

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndex -> unsigned).
   */
  if(_assignmentOption) {
    _shapePositionMap = siteToShapeVertexMap(
      _abstract.permutations.stereopermutations.at(
        _feasible.indices.at(
          _assignmentOption.value()
        )
      ),
      _abstract.canonicalSites
    );
  } else { // Wipe the map
    _shapePositionMap.clear();
  }
}

void AtomStereopermutator::Impl::assignRandom(random::Engine& engine) {
  assign(
    temple::random::pickDiscrete(
      // Map the feasible permutations onto their weights
      temple::map(
        _feasible.indices,
        [&](const unsigned permutationIndex) -> unsigned {
          return _abstract.permutations.weights.at(permutationIndex);
        }
      ),
      engine
    )
  );
}

void AtomStereopermutator::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // RankingInformation changes (lots of atom indices)
  _ranking.applyPermutation(permutation);

  // _centerAtom must change
  _centerAtom = permutation.at(_centerAtom);

  // Neither shape nor assignment change

  /* Although ranking and central atom are implicated in its creation,
   * abstract and feasible's states are independent of atom indices.
   */
}

boost::optional<AtomStereopermutator::PropagatedState> AtomStereopermutator::Impl::propagate(
  const OuterGraph& graph,
  RankingInformation newRanking,
  boost::optional<Symmetry::Shape> shapeOption
) {
  // If nothing changes, nothing changes, and you don't get any internal state
  if(newRanking == _ranking) {
    // TODO no, this isn't right! If the symmetry has changed, then the internal state must change
    return boost::none;
  }

#ifndef NDEBUG
  /* The following is a precondition. It is also rankPriority's postcondition,
   * but propagation may mess with it, so we need to be sure.
   */
  for(const auto& siteAtomList : boost::range::join(_ranking.sites, newRanking.sites)) {
    assert(
      std::is_sorted(
        std::begin(siteAtomList),
        std::end(siteAtomList)
      )
    );
  }
#endif

  /* There are five possible situations that can occur here:
   * - A substituent is removed, and that substituent was sole constituent of a
   *   site, leading to a shape size decrease.
   * - A substituent is removed, but that substituent was part of a site with
   *   multiple constituents. The shape stays the same.
   * - No substituent is removed or added, but there is a ranking change. The
   *   shape stays the same.
   * - A substituent is added, but to an existing site. The shape stays
   *   the same.
   * - A substituent is added that by itself constitutes a new site. The
   *   shape size increases.
   */
  enum class PropagationSituation {
    Unknown,
    SiteRemoval,
    SubstituentRemoval,
    RankingChange,
    SubstituentAddition,
    SiteAddition
  };

  /* Fully determine in which situation we are and all conditional information
   * needed to proceed with each.
   */
  auto situation = PropagationSituation::Unknown;
  int siteCountChange = static_cast<int>(newRanking.sites.size() - _ranking.sites.size());
  boost::optional<unsigned> alteredSiteIndex;
  boost::optional<AtomIndex> alteredSubstituentIndex;

  // Complain if more than one site is changed either way
  if(std::abs(siteCountChange) > 1) {
    throw std::logic_error("propagateGraphChange should only ever handle a single site addition or removal at once");
  }

  auto countSubstituents = [](const RankingInformation::SiteListType& sites) -> unsigned {
    return temple::accumulate(
      sites,
      0u,
      [](const unsigned carry, const auto& site) -> unsigned {
        return carry + site.size();
      }
    );
  };

  int substituentCountChange = static_cast<int>(countSubstituents(newRanking.sites)) - static_cast<int>(countSubstituents(_ranking.sites));
  // Complain if multiple substituents are changed either way
  if(std::abs(substituentCountChange) > 1) {
    throw std::logic_error("propagateGraphChange should only ever handle a single substituent addition or removal at once");
  }

  // Collect all substituents of a ranking and sort them
  auto collectSubstituents = [](const RankingInformation::SiteListType& sites) -> std::vector<AtomIndex> {
    std::vector<AtomIndex> substituents;
    substituents.reserve(2 * sites.size());

    for(const auto& site : sites) {
      for(const AtomIndex& index : site) {
        substituents.push_back(index);
      }
    }

    temple::inplace::sort(substituents);

    return substituents;
  };

  auto determineChangedSubstituentAndSite = [&](
    const RankingInformation& rankingWithMoreSubstituents,
    const RankingInformation& rankingWithFewerSubstituents
  ) {
    // Figure out which substituent has changed, if any
    auto largerSubstituentList = collectSubstituents(rankingWithMoreSubstituents.sites);
    auto smallerSubstituentList = collectSubstituents(rankingWithFewerSubstituents.sites);

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
    alteredSiteIndex = rankingWithMoreSubstituents.getSiteIndexOf(changedSubstituents.front());
  };

  if(substituentCountChange == +1) {
    determineChangedSubstituentAndSite(newRanking, _ranking);

    assert(!newRanking.sites.at(*alteredSiteIndex).empty());
    if(newRanking.sites.at(*alteredSiteIndex).size() > 1) {
      situation = PropagationSituation::SubstituentAddition;
    } else {
      situation = PropagationSituation::SiteAddition;
    }
  } else if(substituentCountChange == -1) {
    determineChangedSubstituentAndSite(_ranking, newRanking);
    assert(!_ranking.sites.at(*alteredSiteIndex).empty());
    if(_ranking.sites.at(*alteredSiteIndex).size() > 1) {
      situation = PropagationSituation::SubstituentRemoval;
    } else {
      situation = PropagationSituation::SiteRemoval;
    }
  } else {
    situation = PropagationSituation::RankingChange;
  }

  /* Decide the new shape */
  Symmetry::Shape newShape = shapeOption.value_or_eval(
    [&]() {
      if(siteCountChange == +1) {
        return up(_shape);
      }

      if(siteCountChange == 0) {
        return _shape;
      }

      assert(siteCountChange == -1);
      assert(alteredSiteIndex);

      /* We can only figure out a mapping if the stereocenter is assigned
       * since otherwise _cache.symmetryPositionMap is empty and the symmetry
       * position being removed is a necessary argument to down.
       */
      if(_assignmentOption) {
        return down(
          _shape,
          _shapePositionMap.at(alteredSiteIndex.value())
        );
      }

      // Just return the most symmetric symmetry of the target size
      return Symmetry::properties::mostSymmetric(
        newRanking.sites.size()
      );
    }
  );

  /* Generate new assignments */
  AbstractStereopermutations newAbstract {
    newRanking,
    newShape
  };

  FeasibleStereopermutations newFeasible {
    newAbstract,
    newShape,
    _centerAtom,
    newRanking,
    graph
  };

  boost::optional<unsigned> newStereopermutationOption;

  /* Now we will attempt to propagate our assignment to the new set of
   * assignments. This is only necessary in the case that the stereopermutator
   * is currently assigned.
   *
   * In the case of an equal number of sites, it is only possible to
   * propagate chiral state if the new number of assignments is smaller or
   * equal to the amount we have currently. This is because say we have AABCDE
   * in octahedral, we carry chiral state, and a ranking change leads to
   * ABCDEF.  There will be multiple ways to split A to F!
   */
  if(_assignmentOption && numStereopermutations() > 1) {
    // A flat map of symmetry position to new site index in the new symmetry
    std::vector<unsigned> sitesAtNewSymmetryPositions;

    /* Site-level changes */
    if(situation == PropagationSituation::SiteAddition) {
      /* Try to get a mapping to the new symmetry. If that returns a Some, try
       * to get a mapping by preservationOption policy. If any of these steps
       * returns boost::none, the whole expression is boost::none.
       */
      auto suitableMappingOption = temple::optionals::flatMap(
        Symmetry::getMapping(
          _shape,
          newShape,
          boost::none
        ),
        [&](const auto& mappingOption) {
          return getIndexMapping(
            mappingOption,
            Options::chiralStatePreservation
          );
        }
      );

      if(suitableMappingOption) {
        /* So now we must transfer the current assignment into the new symmetry
         * and search for it in the set of uniques.
         */
        const auto& symmetryMapping = suitableMappingOption.value();

        // Apply the mapping to get old symmetry positions at their places in the new symmetry
        auto oldSymmetryPositionsInNewSymmetry = Symmetry::properties::applyIndexMapping(
          newShape,
          symmetryMapping
        );

        // Invert symmetryPositionMap to get: site = map.at(symmetryPosition)
        std::vector<unsigned> oldSymmetryPositionToSiteMap(_shapePositionMap.size());
        for(unsigned i = 0 ; i < _shapePositionMap.size(); ++i) {
          oldSymmetryPositionToSiteMap.at(
            _shapePositionMap.at(i)
          ) = i;
        }
        /* We assume the new site is also merely added to the end! This may
         * not always be true
         */
        oldSymmetryPositionToSiteMap.push_back(_shapePositionMap.size());

        // Replace old symmetry positions by their site indices
        sitesAtNewSymmetryPositions = temple::map(
          oldSymmetryPositionsInNewSymmetry,
          [&oldSymmetryPositionToSiteMap](const unsigned oldSymmetryPosition) -> unsigned {
            return oldSymmetryPositionToSiteMap.at(oldSymmetryPosition);
          }
        );
      }
      /* If no mapping can be found that fits to the preservationOption,
       * newStereopermutationOption remains boost::none, and this stereopermutator loses
       * any chiral information it may have had.
       */
    }

    if(situation == PropagationSituation::SiteRemoval) {
      /* Try to get a mapping to the new symmetry. If that returns a Some, try
       * to get a mapping by preservationOption policy. If any of these steps
       * returns boost::none, the whole expression is boost::none.
       */
      auto suitableMappingOptional = temple::optionals::flatMap(
        Symmetry::getMapping(
          _shape,
          newShape,
          /* Last parameter is the deleted symmetry position, which is the
           * symmetry position at which the site being removed is currently at
           */
          _shapePositionMap.at(*alteredSiteIndex)
        ),
        [&](const auto& mappingOptional) {
          return getIndexMapping(mappingOptional, Options::chiralStatePreservation);
        }
      );

      /* It is weird that this block of code works although it is structurally
       * so different from site addition propagation. Is it coincidence?
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

        // Invert symmetryPositionMap to get: site = map.at(symmetryPosition)
        std::vector<unsigned> oldSymmetryPositionToSiteMap(_shapePositionMap.size());
        for(unsigned i = 0 ; i < _shapePositionMap.size(); ++i) {
          oldSymmetryPositionToSiteMap.at(
            _shapePositionMap.at(i)
          ) = i;
        }

        // Transfer indices from current symmetry to new symmetry
        sitesAtNewSymmetryPositions.resize(Symmetry::size(newShape));
        for(unsigned i = 0; i < Symmetry::size(newShape); ++i) {
          // i is a symmetry position
          sitesAtNewSymmetryPositions.at(i) = oldSymmetryPositionToSiteMap.at(
            symmetryMapping.at(i)
          );
        }
      }
    }

    /* Substituent-level changes */
    if(situation == PropagationSituation::SubstituentAddition) {
      /* Sort sites' constituting atom indices in both rankings so we can use
       * lexicographical comparison to create a mapping. Add the
       * newSubstituentIndex to the old ranking's appropriate site so the
       * mapping is 1:1.
       */
      temple::TinySet<AtomIndex>::checked_insert(
        _ranking.sites.at(*alteredSiteIndex),
        *alteredSubstituentIndex
      );

      // Generate a flat map of old site indices to new site indices
      auto siteMapping = temple::map(
        _ranking.sites,
        [&newRanking](const auto& siteAtomList) -> unsigned {
          auto findIter = std::find(
            newRanking.sites.begin(),
            newRanking.sites.end(),
            siteAtomList
          );

          assert(findIter != newRanking.sites.end());

          return findIter - newRanking.sites.begin();
        }
      );

      // Write the found mapping into sitesAtNewSymmetryPositions
      sitesAtNewSymmetryPositions.resize(Symmetry::size(newShape));
      for(unsigned i = 0; i < siteMapping.size(); ++i) {
        sitesAtNewSymmetryPositions.at(i) = siteMapping.at(
          _shapePositionMap.at(i)
        );
      }

      /* Revert the change to _ranking from earlier for easy mapping so that
       * the returned prior state is unchanged.
       */
      temple::inplace::remove(
        _ranking.sites.at(*alteredSiteIndex),
        *alteredSubstituentIndex
      );
    }

    if(situation == PropagationSituation::SubstituentRemoval) {
      /* Sort sites in the old ranking and new so we can use lexicographical
       * comparison to figure out a mapping
       */
      temple::inplace::remove(
        _ranking.sites.at(*alteredSiteIndex),
        *alteredSubstituentIndex
      );

      // Calculate the mapping from old sites to new ones
      auto siteMapping = temple::map(
        _ranking.sites,
        [&newRanking](const auto& siteAtomList) -> unsigned {
          auto findIter = std::find(
            newRanking.sites.begin(),
            newRanking.sites.end(),
            siteAtomList
          );

          assert(findIter != newRanking.sites.end());

          return findIter - newRanking.sites.begin();
        }
      );

      sitesAtNewSymmetryPositions.resize(Symmetry::size(newShape));
      // Transfer sites to new mapping
      for(unsigned i = 0; i < siteMapping.size(); ++i) {
        sitesAtNewSymmetryPositions.at(i) = siteMapping.at(
          _shapePositionMap.at(i)
        );
      }

      /* Revert the change to _ranking from earlier for easy mapping so that
       * the returned prior state is unchanged.
       */
      temple::TinySet<AtomIndex>::checked_insert(
        _ranking.sites.at(*alteredSiteIndex),
        *alteredSubstituentIndex
      );
    }

    /* Ranking-level change */
    if(
      situation == PropagationSituation::RankingChange
      && (
        newAbstract.permutations.stereopermutations.size()
        <= _abstract.permutations.stereopermutations.size()
      )
    ) {
      const auto& currentStereopermutation = _abstract.permutations.stereopermutations.at(
        _feasible.indices.at(
          _assignmentOption.value()
        )
      );

      // Replace the characters by their corresponding indices from the old ranking
      sitesAtNewSymmetryPositions = shapeVertexToSiteIndexMap(
        currentStereopermutation,
        _abstract.canonicalSites
      );
    }

    if(!sitesAtNewSymmetryPositions.empty()) {
      // Replace the site indices by their new ranking characters
      auto newStereopermutationCharacters = AbstractStereopermutations::makeStereopermutationCharacters(
        newAbstract.canonicalSites,
        newAbstract.symbolicCharacters,
        sitesAtNewSymmetryPositions
      );

      // Create a new assignment with those characters
      auto trialStereopermutation = stereopermutation::Stereopermutation(
        newStereopermutationCharacters,
        newAbstract.selfReferentialLinks
      );

      // Generate all rotations of this trial assignment
      auto allTrialRotations = trialStereopermutation.generateAllRotations(newShape);

      // Find out which of the new assignments has a rotational equivalent
      for(unsigned i = 0; i < newAbstract.permutations.stereopermutations.size(); ++i) {
        if(allTrialRotations.count(newAbstract.permutations.stereopermutations.at(i)) > 0) {
          newStereopermutationOption = i;
          break;
        }
      }
    }
  }

  /* If we have discovered a stereopermutation to map to, we must still figure
   * out its assignment index (assuming it is feasible)
   */
  auto newAssignmentOption = temple::optionals::flatMap(
    newStereopermutationOption,
    [&](const unsigned stereopermutationIndex) -> boost::optional<unsigned> {
      auto assignmentFindIter = std::find(
        std::begin(newFeasible.indices),
        std::end(newFeasible.indices),
        stereopermutationIndex
      );

      if(assignmentFindIter == std::end(newFeasible.indices)) {
        // The stereopermutation is infeasible
        return boost::none;
      }

      return assignmentFindIter - std::begin(newFeasible.indices);
    }
  );

  // Extract old state from the class
  auto oldStateTuple = std::make_tuple(
    std::move(_ranking),
    std::move(_abstract),
    std::move(_feasible),
    std::move(_assignmentOption)
  );

  // Overwrite the class state
  _shape = newShape;
  _ranking = std::move(newRanking);
  _abstract = std::move(newAbstract);
  _feasible = std::move(newFeasible);
  assign(newAssignmentOption);

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
  for(auto& equalPrioritySet : _ranking.substituentRanking) {
    for(auto& index : equalPrioritySet) {
      updateIndexInplace(index);
    }
  }

  for(auto& siteAtomList : _ranking.sites) {
    for(auto& atomIndex : siteAtomList) {
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

const AbstractStereopermutations& AtomStereopermutator::Impl::getAbstract() const {
  return _abstract;
}

const FeasibleStereopermutations& AtomStereopermutator::Impl::getFeasible() const {
  return _feasible;
}

const RankingInformation& AtomStereopermutator::Impl::getRanking() const {
  return _ranking;
}

Symmetry::Shape AtomStereopermutator::Impl::getShape() const {
  return _shape;
}

const std::vector<unsigned>& AtomStereopermutator::Impl::getShapePositionMap() const {
  if(_assignmentOption == boost::none) {
    throw std::logic_error(
      "The AtomStereopermutator is unassigned, sites are not assigned to "
      "symmetry positions"
    );
  }

  return _shapePositionMap;
}

void AtomStereopermutator::Impl::fit(
  const OuterGraph& graph,
  const AngstromWrapper& angstromWrapper
) {
  const unsigned S = Symmetry::size(_shape);
  assert(S == _ranking.sites.size());

  // For all atoms making up a site, decide on the spatial average position
  Eigen::Matrix<double, 3, Eigen::Dynamic> sitePositions(3, S);
  for(unsigned i = 0; i < S; ++i) {
    sitePositions.col(i) = cartesian::averagePosition(angstromWrapper.positions, _ranking.sites.at(i));
  }

  std::vector<Symmetry::Shape> excludeShapes;

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
          cartesian::angle(
            sitePositions.col(i),
            angstromWrapper.positions.row(_centerAtom),
            sitePositions.col(j)
          )
        );
      }
    }

    std::sort(std::begin(angles), std::end(angles));

    const double tau = Symmetry::tau(angles);

    if(S == 4) {
      /* Thresholds
       * - τ₄' = 0 -> Symmetry is square planar
       * - τ₄' = 0.24 -> Symmetry is seesaw
       * - τ₄' = 1 -> Symmetry is tetrahedral
       */
      if(tau < 0.12) {
        // Symmetry is square planar
        excludeShapes.push_back(Symmetry::Shape::Disphenoid);
        excludeShapes.push_back(Symmetry::Shape::Tetrahedron);
      } else if(0.12 <= tau && tau < 0.62) {
        excludeShapes.push_back(Symmetry::Shape::Square);
        // Symmetry is seesaw
        excludeShapes.push_back(Symmetry::Shape::Tetrahedron);
      } else if(0.62 <= tau) {
        excludeShapes.push_back(Symmetry::Shape::Square);
        excludeShapes.push_back(Symmetry::Shape::Disphenoid);
        // Symmetry is tetrahedral
      }
    } else if(S == 5) {
      /* Thresholds:
       * - τ₅ = 0 -> Symmetry is square pyramidal
       * - τ₅ = 1 -> Symmetry is trigonal bipyramidal
       */

      if(tau < 0.5) {
        excludeShapes.push_back(Symmetry::Shape::TrigonalBipyramid);
      } else if(tau > 0.5) {
        excludeShapes.push_back(Symmetry::Shape::SquarePyramid);
      }
    }
  }

  // Save stereopermutator state to return to if no fit is viable
  const Symmetry::Shape priorShape = _shape;
  const boost::optional<unsigned> priorStereopermutation  = _assignmentOption;

  const Symmetry::Shape initialShape {Symmetry::Shape::Line};
  const unsigned initialStereopermutation = 0;
  const double initialPenalty = 100;

  Symmetry::Shape bestShape = initialShape;
  unsigned bestStereopermutation = initialStereopermutation;
  double bestPenalty = initialPenalty;
  unsigned bestStereopermutationMultiplicity = 1;

  auto excludesContains = temple::makeContainsPredicate(excludeShapes);

  // Cycle through all symmetries
  for(const auto& shape : Symmetry::allShapes) {
    // Skip any shapes of different size
    if(Symmetry::size(shape) != S || excludesContains(shape)) {
      continue;
    }

    // Change the shape of the AtomStereopermutator
    setShape(shape, graph);

    const unsigned assignmentCount = numAssignments();

    for(unsigned assignment = 0; assignment < assignmentCount; ++assignment) {
      // Assign the stereopermutator
      assign(assignment);

      const double angleDeviations = temple::sum(
        temple::adaptors::transform(
          temple::adaptors::allPairs(
            temple::adaptors::range(Symmetry::size(_shape))
          ),
          [&](const unsigned siteI, const unsigned siteJ) -> double {
            return std::fabs(
              cartesian::angle(
                sitePositions.col(siteI),
                angstromWrapper.positions.row(_centerAtom),
                sitePositions.col(siteJ)
              ) - angle(siteI, siteJ)
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
            temple::adaptors::range(Symmetry::size(_shape))
          ),
          [&](const unsigned siteI, const unsigned siteJ) -> double {
            return std::fabs(
              // siteI - siteJ 1-3 distance from positions
              cartesian::distance(
                sitePositions.col(siteI),
                sitePositions.col(siteJ)
              )
              // idealized 1-3 distance from
              - CommonTrig::lawOfCosines(
                // i-j 1-2 distance from positions
                cartesian::distance(
                  sitePositions.col(siteI),
                  angstromWrapper.positions.row(_centerAtom)
                ),
                // j-k 1-2 distance from positions
                cartesian::distance(
                  angstromWrapper.positions.row(_centerAtom),
                  sitePositions.col(siteJ)
                ),
                // idealized Stereopermutator angle
                angle(siteI, siteJ)
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
          minimalChiralConstraints(),
          [&](const auto& minimalPrototype) -> double {
            auto fetchPosition = [&](const boost::optional<unsigned>& siteIndexOptional) -> Eigen::Vector3d {
              if(siteIndexOptional) {
                return sitePositions.col(siteIndexOptional.value());
              }

              return angstromWrapper.positions.row(_centerAtom);
            };

            double volume = cartesian::adjustedSignedVolume(
              fetchPosition(minimalPrototype[0]),
              fetchPosition(minimalPrototype[1]),
              fetchPosition(minimalPrototype[2]),
              fetchPosition(minimalPrototype[3])
            );

            // minimalChiralConstraints() supplies only Positive targets
            if(volume < 0) {
              return 1;
            }

            return 0;
          }
        )
      );

      double fitPenalty = (
        angleDeviations
        + oneThreeDistanceDeviations
        + chiralityDeviations
      );


#ifndef NDEBUG
      Log::log(Log::Particulars::AtomStereopermutatorFit)
        << Symmetry::nameIndex(shape)
        << ", " << assignment
        << ", " << std::setprecision(4) << std::fixed
        << angleDeviations << ", "
        << oneThreeDistanceDeviations << ", "
        << chiralityDeviations
        << std::endl;
#endif

      if(fitPenalty < bestPenalty) {
        bestShape = shape;
        bestStereopermutation = assignment;
        bestPenalty = fitPenalty;
        bestStereopermutationMultiplicity = 1;
      } else if(fitPenalty == bestPenalty) {
        // Assume that IF we have multiplicity, it's from the same symmetry
        assert(bestShape == shape);
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
    bestShape == initialShape
    && bestStereopermutation == initialStereopermutation
    && bestPenalty == initialPenalty
  ) {
    // Return to prior
    setShape(priorShape, graph);
    assign(priorStereopermutation);
  } else {
    // Set to best fit
    setShape(bestShape, graph);

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
  assert(!_shapePositionMap.empty());

  return Symmetry::angleFunction(_shape)(
    _shapePositionMap.at(i),
    _shapePositionMap.at(j)
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
    return _feasible.indices.at(_assignmentOption.value());
  }

  return boost::none;
}

std::vector<AtomStereopermutator::MinimalChiralConstraint>
AtomStereopermutator::Impl::minimalChiralConstraints(bool enforce) const {
  std::vector<
    std::array<boost::optional<unsigned>, 4>
  > precursors;

  /* It only makes sense to emit these minimal representations of chiral
   * constraints if the stereopermutator is assigned.
   *
   * As long as the permutator is achiral, there is little reason to burden
   * refinement with error terms that achieve nothing for the final
   * conformation. Yet, if any sort of BondStereopermutator is placed on even
   * an achiral AtomStereopermutator, assumptions are made about the spatial
   * placement of symmetry positions.
   *
   * Since refinement has many terms, it can be difficult to reliably minimize
   * these assumptions on symmetry position spatial locations against partially
   * relaxed distance constraints without the clear ordering priority of chiral
   * constraints enforcing your spatial reasoning assumptions.
   *
   * So for now, until a better solution arises, emit chiral constraints if a
   * bond stereopermutator relies upon symmetry position assumptions even if
   * this stereopermutator might be achiral.
   */
  if(
    _assignmentOption
    && (
      numStereopermutations() > 1
      || enforce
    )
  ) {

    /* Invert _neighborSymmetryPositionMap, we need a mapping of
     *  (position in symmetry) -> atom index
     */
    auto symmetryPositionToSiteIndexMap = shapeVertexToSiteIndexMap(
      _abstract.permutations.stereopermutations.at(
        _feasible.indices.at(
          _assignmentOption.value()
        )
      ),
      _abstract.canonicalSites
    );

    // Get list of tetrahedra from symmetry
    const auto& tetrahedraList = Symmetry::tetrahedra(_shape);

    precursors.reserve(tetrahedraList.size());
    for(const auto& tetrahedron : tetrahedraList) {
      /* Replace indices (represent positions within the symmetry) with the
       * site index at that position from the inverted map
       */

      // Make a minimal sequence from it
      precursors.push_back(
        temple::map(
          tetrahedron,
          [&](const boost::optional<unsigned>& indexOptional) -> boost::optional<unsigned> {
            if(indexOptional) {
              return symmetryPositionToSiteIndexMap.at(
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
    + std::to_string(_centerAtom) + " ("s + Symmetry::name(_shape) +", "s;

  const auto& characters = _abstract.symbolicCharacters;
  std::copy(
    characters.begin(),
    characters.end(),
    std::back_inserter(returnString)
  );

  for(const auto& link : _abstract.selfReferentialLinks) {
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
    "A-"s + std::to_string(static_cast<unsigned>(_shape))
    + "-"s + std::to_string(numStereopermutations())
    + "-"s + (
      indexOfPermutation()
      ? std::to_string(indexOfPermutation().value())
      : "u"s
    )
  );
}

unsigned AtomStereopermutator::Impl::numAssignments() const {
  return _feasible.indices.size();
}

unsigned AtomStereopermutator::Impl::numStereopermutations() const {
  return _abstract.permutations.stereopermutations.size();
}

void AtomStereopermutator::Impl::setShape(
  const Symmetry::Shape shape,
  const OuterGraph& graph
) {
  if(_shape == shape) {
    // If the symmetry doesn't actually change, then nothing does
    return;
  }

  _shape = shape;

  _abstract = AbstractStereopermutations {
    _ranking,
    _shape
  };

  _feasible = FeasibleStereopermutations {
    _abstract,
    _shape,
    _centerAtom,
    _ranking,
    graph
  };

  // Dis-assign the stereopermutator
  assign(boost::none);
}

} // namespace molassembler

} // namespace Scine
