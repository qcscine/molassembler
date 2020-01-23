/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/AtomStereopermutatorImpl.h"

#include "boost/range/join.hpp"
#include "shapes/Properties.h"
#include "shapes/PropertyCaching.h"
#include "shapes/ContinuousMeasures.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "stereopermutation/Manipulation.h"
#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"
#include "temple/Functor.h"
#include "temple/Optionals.h"
#include "temple/OrderedPair.h"
#include "temple/Random.h"
#include "temple/TinySet.h"
#include "temple/constexpr/Math.h"
#include "temple/constexpr/Numeric.h"

#include "molassembler/Cycles.h"
#include "molassembler/Detail/BuildTypeSwitch.h"
#include "molassembler/Detail/Cartesian.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/Graph/PrivateGraph.h"
#include "molassembler/Log.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/Stereopermutators/ShapeVertexMaps.h"

#include "Utils/Geometry/ElementInfo.h"

namespace Scine {

namespace molassembler {

/* Static functions */

namespace detail {

shapes::Shape pickTransition(
  const shapes::Shape shape,
  const unsigned T,
  boost::optional<unsigned> removedSymmetryPositionOptional
) {
  boost::optional<shapes::properties::SymmetryTransitionGroup> bestTransition;
  std::vector<shapes::Shape> propositions;

  auto replaceOrAdd = [&](
    shapes::Shape newName,
    const shapes::properties::SymmetryTransitionGroup& transitionGroup
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
  for(const shapes::Shape propositionalShape : shapes::allShapes) {
    if(shapes::size(propositionalShape) != T) {
      continue;
    }

    if(auto transitionOptional = shapes::getMapping(shape, propositionalShape, removedSymmetryPositionOptional)) {
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
    return shapes::properties::mostSymmetric(T);
  }

  return shapes::properties::mostSymmetric(std::move(propositions));
}

std::pair<shapes::Shape, std::vector<unsigned>> classifyShape(const Eigen::Matrix<double, 3, Eigen::Dynamic>& sitePositions) {
  const unsigned S = sitePositions.cols() - 1;
  auto normalized = shapes::continuous::normalize(sitePositions);

  std::vector<shapes::Shape> viableShapes;
  for(const shapes::Shape shape : shapes::allShapes) {
    if(shapes::size(shape) == S) {
      viableShapes.push_back(shape);
    }
  }
  const unsigned shapesCount = viableShapes.size();
  std::vector<shapes::continuous::ShapeResult> shapeMeasureResults (shapesCount);

#pragma omp parallel for
  for(unsigned i = 0; i < shapesCount; ++i) {
    const shapes::Shape candidateShape = viableShapes[i];
    shapeMeasureResults[i] = shapes::continuous::shapeCentroidLast(normalized, candidateShape);

    // Bias shape classification against some shapes
    if(candidateShape == shapes::Shape::TrigonalPyramid) {
      shapeMeasureResults[i].measure *= 4;
    } else if(candidateShape == shapes::Shape::Seesaw) {
      shapeMeasureResults[i].measure *= 2;
    }
  }

  // Ensure centroids are mapped against one another
  assert(
    temple::all_of(
      shapeMeasureResults,
      [](const auto& shapeResult) -> bool {
        assert(!shapeResult.mapping.empty());
        return shapeResult.mapping.back() == shapeResult.mapping.size() - 1;
      }
    )
  );

  const auto minElementIter = std::min_element(
    std::begin(shapeMeasureResults),
    std::end(shapeMeasureResults),
    [](const auto& a, const auto& b) -> bool {
      return a.measure < b.measure;
    }
  );

  const unsigned minimalShapeIndex = minElementIter - std::begin(shapeMeasureResults);
  const shapes::Shape minimalShape = viableShapes.at(minimalShapeIndex);
  return std::make_pair(
    minimalShape,
    std::move(minElementIter->mapping)
  );
}

} // namespace detail

shapes::Shape AtomStereopermutator::Impl::up(const shapes::Shape shape) {
  return detail::pickTransition(shape, shapes::size(shape) + 1, boost::none);
}

shapes::Shape AtomStereopermutator::Impl::down(const shapes::Shape shape, const unsigned removedShapePosition) {
  return detail::pickTransition(shape, shapes::size(shape) - 1, removedShapePosition);
}

boost::optional<std::vector<unsigned>> AtomStereopermutator::Impl::getIndexMapping(
  const shapes::properties::SymmetryTransitionGroup& mappingsGroup,
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

bool AtomStereopermutator::Impl::thermalized(
  const Graph& graph,
  const AtomIndex centerAtom,
  const shapes::Shape shape,
  const RankingInformation& ranking,
  const TemperatureRegime temperature
) {
  if(temperature == TemperatureRegime::Low) {
    return false;
  }

  /* Nitrogen atom inversion */
  constexpr unsigned nitrogenZ = 7;
  bool isNitrogenIsotope = Utils::ElementInfo::Z(graph.elementType(centerAtom)) == nitrogenZ;

  if(
    isNitrogenIsotope
    && shape == shapes::Shape::VacantTetrahedron
  ) {
    // Generally thermalized, except if in a small cycle
    if(
      temple::any_of(ranking.links,
        [](const LinkInformation& link) {
          return link.cycleSequence.size() <= 4;
        }
      )
    ) {
      return false;
    }

    return true;
  }


  // Berry pseudorotation and Bartell mechanism
  if(
    ranking.links.empty() && (
      shape == shapes::Shape::PentagonalBipyramid
      || shape == shapes::Shape::TrigonalBipyramid
    )
  ) {
    return true;
  }

  return false;
}

/* Constructors */
AtomStereopermutator::Impl::Impl(
  const Graph& graph,
  // The symmetry of this Stereopermutator
  const shapes::Shape shape,
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
    _shapePositionMap {},
    _thermalized {thermalized(
      graph,
      centerAtom,
      shape,
      ranking,
      Options::temperatureRegime
    )}
{}

/* Modification */
void AtomStereopermutator::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment && assignment.value() >= _feasible.indices.size()) {
    throw std::out_of_range("Supplied assignment index is out of range");
  }

  // Store new assignment
  _assignmentOption = std::move(assignment);

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndex -> unsigned).
   */
  if(_assignmentOption) {
    _shapePositionMap = siteToShapeVertexMap(
      _abstract.permutations.list.at(
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
  const unsigned A = numAssignments();
  if(A == 0) {
    throw std::logic_error("Cannot randomly assign a stereopermutator without feasible stereopermutations");
  }

  if(A == 1) {
    assign(0);
  } else {
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
}

void AtomStereopermutator::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // RankingInformation changes (lots of atom indices)
  _ranking.applyPermutation(permutation);

  // _centerAtom must change
  _centerAtom = permutation.at(_centerAtom);

  /* Neither shape nor assignment change. Also, although ranking and central
   * atom are implicated in their creation, abstract and feasible's states are
   * independent of atom indices.
   */
}

boost::optional<AtomStereopermutator::PropagatedState> AtomStereopermutator::Impl::propagate(
  const Graph& graph,
  RankingInformation newRanking,
  boost::optional<shapes::Shape> shapeOption
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
  shapes::Shape newShape = shapeOption.value_or_eval(
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
      return shapes::properties::mostSymmetric(
        newRanking.sites.size()
      );
    }
  );

  /* Generate new assignments */
  stereopermutators::Abstract newAbstract {
    newRanking,
    newShape
  };

  stereopermutators::Feasible newFeasible {
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
    std::vector<unsigned> sitesAtNewShapeVertices;

    /* Site-level changes */
    if(situation == PropagationSituation::SiteAddition) {
      /* Try to get a mapping to the new symmetry. If that returns a Some, try
       * to get a mapping by preservationOption policy. If any of these steps
       * returns boost::none, the whole expression is boost::none.
       */
      auto suitableMappingOption = temple::optionals::flatMap(
        shapes::getMapping(_shape, newShape, boost::none),
        [](const auto& mappingOption) {
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
        auto oldSymmetryPositionsInNewSymmetry = shapes::properties::applyIndexMapping(
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
        sitesAtNewShapeVertices = temple::map(
          oldSymmetryPositionsInNewSymmetry,
          temple::functor::at(oldSymmetryPositionToSiteMap)
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
        shapes::getMapping(
          _shape,
          newShape,
          /* Last parameter is the deleted symmetry position, which is the
           * symmetry position at which the site being removed is currently at
           */
          _shapePositionMap.at(*alteredSiteIndex)
        ),
        [](const auto& mappingOptional) {
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

        // Invert shapePositionMap to get: site = map.at(shapeVertex)
        std::vector<unsigned> oldShapeVertexToSiteMap(_shapePositionMap.size());
        for(unsigned i = 0 ; i < _shapePositionMap.size(); ++i) {
          oldShapeVertexToSiteMap.at(
            _shapePositionMap.at(i)
          ) = i;
        }

        // Transfer site indices from current shape to new shape
        sitesAtNewShapeVertices.resize(shapes::size(newShape));
        for(unsigned i = 0; i < shapes::size(newShape); ++i) {
          // i is a symmetry position
          sitesAtNewShapeVertices.at(i) = oldShapeVertexToSiteMap.at(
            symmetryMapping.at(i)
          );
        }
      }
    }

    /* Substituent-level changes */
    if(situation == PropagationSituation::SubstituentAddition) {
      /* Add the newSubstituentIndex to the old ranking's appropriate site so
       * the mapping is 1:1.
       */
      RankingInformation::SiteListType oldSites = _ranking.sites;
      temple::TinySet<AtomIndex>::checked_insert(
        oldSites.at(*alteredSiteIndex),
        *alteredSubstituentIndex
      );

      // Generate a flat map of old site indices to new site indices
      auto siteMapping = temple::map(oldSites, temple::functor::indexIn(newRanking.sites));

      // Write the found mapping into sitesAtNewShapeVertices
      sitesAtNewShapeVertices.resize(shapes::size(newShape));
      for(unsigned i = 0; i < siteMapping.size(); ++i) {
        sitesAtNewShapeVertices.at(i) = siteMapping.at(
          _shapePositionMap.at(i)
        );
      }
    }

    if(situation == PropagationSituation::SubstituentRemoval) {
      /* Sort sites in the old ranking and new so we can use lexicographical
       * comparison to figure out a mapping
       */
      RankingInformation::SiteListType oldSites = _ranking.sites;
      temple::inplace::remove(
        oldSites.at(*alteredSiteIndex),
        *alteredSubstituentIndex
      );

      // Calculate the mapping from old sites to new ones
      auto siteMapping = temple::map(oldSites, temple::functor::indexIn(newRanking.sites));

      sitesAtNewShapeVertices.resize(shapes::size(newShape));
      // Transfer sites to new mapping
      for(unsigned i = 0; i < siteMapping.size(); ++i) {
        sitesAtNewShapeVertices.at(i) = siteMapping.at(
          _shapePositionMap.at(i)
        );
      }
    }

    /* Ranking-level change */
    if(
      situation == PropagationSituation::RankingChange
      && (
        newAbstract.permutations.list.size()
        <= _abstract.permutations.list.size()
      )
    ) {
      const auto& currentStereopermutation = _abstract.permutations.list.at(
        _feasible.indices.at(
          _assignmentOption.value()
        )
      );

      // Replace the characters by their corresponding indices from the old ranking
      sitesAtNewShapeVertices = shapeVertexToSiteIndexMap(
        currentStereopermutation,
        _abstract.canonicalSites
      );
    }

    if(!sitesAtNewShapeVertices.empty()) {
      // Replace the site indices by their new ranking characters
      auto newStereopermutationCharacters = stereopermutators::Abstract::makeStereopermutationCharacters(
        newAbstract.canonicalSites,
        newAbstract.symbolicCharacters,
        sitesAtNewShapeVertices
      );

      // Create a new assignment with those characters
      auto trialStereopermutation = stereopermutation::Stereopermutation(
        newStereopermutationCharacters,
        newAbstract.selfReferentialLinks
      );

      // Generate all rotations of this trial assignment
      auto allTrialRotations = stereopermutation::generateAllRotations(trialStereopermutation, newShape);

      // Find out which of the new assignments has a rotational equivalent
      for(unsigned i = 0; i < newAbstract.permutations.list.size(); ++i) {
        auto findIter = std::find(
          std::begin(allTrialRotations),
          std::end(allTrialRotations),
          newAbstract.permutations.list.at(i)
        );
        if(findIter != std::end(allTrialRotations)) {
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
  _thermalized = thermalized(
    graph,
    _centerAtom,
    _shape,
    _ranking,
    Options::temperatureRegime
  );
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
      index = PrivateGraph::removalPlaceholder;
    }
  };

  auto updateIndex = [&removedIndex](const AtomIndex index) -> AtomIndex {
    if(index > removedIndex) {
      return index - 1;
    }

    if(index == removedIndex) {
      return PrivateGraph::removalPlaceholder;
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

const stereopermutators::Abstract& AtomStereopermutator::Impl::getAbstract() const {
  return _abstract;
}

const stereopermutators::Feasible& AtomStereopermutator::Impl::getFeasible() const {
  return _feasible;
}

const RankingInformation& AtomStereopermutator::Impl::getRanking() const {
  return _ranking;
}

shapes::Shape AtomStereopermutator::Impl::getShape() const {
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
  const Graph& graph,
  const AngstromWrapper& angstromWrapper
) {
  const unsigned S = shapes::size(_shape);
  assert(S == _ranking.sites.size());

  // Save stereopermutator state to return to if no fit is viable
  const shapes::Shape priorShape = _shape;
  const boost::optional<unsigned> priorStereopermutation  = _assignmentOption;

  // For all atoms making up a site, decide on the spatial average position
  Eigen::Matrix<double, 3, Eigen::Dynamic> sitePositions(3, S + 1);
  for(unsigned i = 0; i < S; ++i) {
    sitePositions.col(i) = cartesian::averagePosition(angstromWrapper.positions, _ranking.sites.at(i));
  }
  // Add the putative center
  sitePositions.col(S) = angstromWrapper.positions.row(_centerAtom);

  // Classify the shape and set it
  shapes::Shape fittedShape;
  std::vector<unsigned> matchingMapping;
  std::tie(fittedShape, matchingMapping) = detail::classifyShape(sitePositions);
  /* Drop the centroid, making the remaining sequence viable as an index mapping
   * within the set of shape rotations, i.e. for use in stereopermutation
   * functionality
   */
  matchingMapping.pop_back();

  setShape(fittedShape, graph);

  /* Ok, so: We have a mapping from site positions to shape vertices from the
   * shape classification algorithm. Site positions are ordered identically to
   * the sites, so we essentially have a mapping from sites to shape vertices.
   *
   * Now we need to find the (feasible) stereopermutation that matches it.
   *
   * Can we transform the site -> shape vertex mapping into a
   * stereopermutation, generate all its rotations and then set-membership
   * check each feasible permutation?
   */
  auto soughtStereopermutation = stereopermutationFromSiteToShapeVertexMap(
    matchingMapping,
    _ranking.links,
    _abstract.canonicalSites
  );

  auto soughtRotations = stereopermutation::generateAllRotations(soughtStereopermutation, fittedShape);

  /* Although the stereopermutations in _abstract are sorted by their index of
   * permutation and we could potentially binary search all of our sought
   * rotations within that set, this linear search will never be time-critical
   * as long as we use continuous shape measure-based shape classification.
   */
  boost::optional<unsigned> foundStereopermutation;
  const unsigned A = _feasible.indices.size();
  for(unsigned a = 0; a < A; ++a) {
    const auto& feasiblePermutation = _abstract.permutations.list.at(
      _feasible.indices.at(a)
    );
    auto findIter = std::find(
      std::begin(soughtRotations),
      std::end(soughtRotations),
      feasiblePermutation
    );
    if(findIter != std::end(soughtRotations)) {
      foundStereopermutation = a;
      break;
    }
  }

  /* In case NO assignments could be found, return to the prior state.
   *
   * This guards against situations in which predicates in
   * uniques could lead no assignments to be returned, such as
   * in e.g. square-planar AAAB with {0, 3}, {1, 3}, {2, 3} with removal of
   * trans-spanning groups. In that situation, all possible assignments are
   * trans-spanning and uniques is an empty vector.
   *
   * At the moment, this predicate is disabled, so no such issues should arise.
   * Just being safe.
   */
  if(foundStereopermutation == boost::none) {
    setShape(priorShape, graph);
    assign(priorStereopermutation);
  } else {
    assign(*foundStereopermutation);
  }
}

/* Information */
double AtomStereopermutator::Impl::angle(
  const unsigned i,
  const unsigned j
) const {
  assert(i != j);
  assert(!_shapePositionMap.empty());

  return shapes::angleFunction(_shape)(
    _shapePositionMap.at(i),
    _shapePositionMap.at(j)
  );
}

boost::optional<unsigned> AtomStereopermutator::Impl::assigned() const {
  if(_thermalized) {
    return temple::optionals::map(_assignmentOption, [](unsigned /* a */) { return 0u; });
  }

  return _assignmentOption;
}

AtomIndex AtomStereopermutator::Impl::centralIndex() const {
  return _centerAtom;
}

boost::optional<unsigned> AtomStereopermutator::Impl::indexOfPermutation() const {
  if(_thermalized) {
    return temple::optionals::map(_assignmentOption, [](unsigned /* a */) { return 0u; });
  }

  return temple::optionals::map(
    _assignmentOption,
    [&](unsigned a) { return _feasible.indices.at(a); }
  );
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
    && (numStereopermutations() > 1 || enforce)
  ) {
    /* Invert _neighborSymmetryPositionMap, we need a mapping of
     *  (position in symmetry) -> atom index
     */
    auto symmetryPositionToSiteIndexMap = shapeVertexToSiteIndexMap(
      _abstract.permutations.list.at(
        _feasible.indices.at(
          _assignmentOption.value()
        )
      ),
      _abstract.canonicalSites
    );

    // Get list of tetrahedra from symmetry
    const auto& tetrahedraList = shapes::tetrahedra(_shape);

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
    + std::to_string(_centerAtom) + " ("s + shapes::name(_shape) +", "s;

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
  if(_thermalized) {
    return 1;
  }

  return _feasible.indices.size();
}

unsigned AtomStereopermutator::Impl::numStereopermutations() const {
  if(_thermalized) {
    return 1;
  }

  return _abstract.permutations.list.size();
}

void AtomStereopermutator::Impl::setShape(
  const shapes::Shape shape,
  const Graph& graph
) {
  if(_shape == shape) {
    // If the symmetry doesn't actually change, then nothing does
    return;
  }

  _shape = shape;

  _abstract = stereopermutators::Abstract {
    _ranking,
    _shape
  };

  _feasible = stereopermutators::Feasible {
    _abstract,
    _shape,
    _centerAtom,
    _ranking,
    graph
  };

  _thermalized = thermalized(
    graph,
    _centerAtom,
    _shape,
    _ranking,
    Options::temperatureRegime
  );

  // Dis-assign the stereopermutator
  assign(boost::none);
}

} // namespace molassembler

} // namespace Scine
