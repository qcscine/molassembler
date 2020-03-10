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
#include "molassembler/Stereopermutators/RankingMapping.h"

#include "Utils/Geometry/ElementInfo.h"

using namespace std::string_literals;
using namespace std::placeholders;

namespace Scine {
namespace molassembler {

/* Static functions */

namespace {

shapes::Shape pickTransition(
  const shapes::Shape shape,
  const unsigned T,
  boost::optional<shapes::Vertex> removedVertexOptional
) {
  boost::optional<shapes::properties::ShapeTransitionGroup> bestTransition;
  std::vector<shapes::Shape> propositions;

  auto replaceOrAdd = [&](
    shapes::Shape newName,
    const shapes::properties::ShapeTransitionGroup& transitionGroup
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

    if(auto transitionOptional = shapes::getMapping(shape, propositionalShape, removedVertexOptional)) {
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

std::pair<shapes::Shape, std::vector<shapes::Vertex>> classifyShape(const Eigen::Matrix<double, 3, Eigen::Dynamic>& sitePositions) {
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

} // namespace

shapes::Shape AtomStereopermutator::Impl::up(const shapes::Shape shape) {
  return pickTransition(shape, shapes::size(shape) + 1, boost::none);
}

shapes::Shape AtomStereopermutator::Impl::down(const shapes::Shape shape, const shapes::Vertex removedVertex) {
  return pickTransition(shape, shapes::size(shape) - 1, removedVertex);
}

boost::optional<std::vector<shapes::Vertex>> AtomStereopermutator::Impl::selectTransitionMapping(
  const shapes::properties::ShapeTransitionGroup& mappingsGroup,
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
  constexpr unsigned nitrogenZ = Utils::ElementInfo::Z(Utils::ElementType::N);
  bool isNitrogenIsotope = Utils::ElementInfo::Z(graph.elementType(centerAtom)) == nitrogenZ;

  if(
    isNitrogenIsotope
    && shape == shapes::Shape::VacantTetrahedron
  ) {
    // Generally thermalized, except if in a small cycle
    if(
      temple::any_of(
        ranking.links,
        [](const RankingInformation::Link& link) {
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
) : centerAtom_ {centerAtom},
    shape_ {shape},
    ranking_ {std::move(ranking)},
    abstract_ {ranking_, shape_},
    feasible_ {abstract_, shape_, centerAtom_, ranking_, graph},
    assignmentOption_ {boost::none},
    shapePositionMap_ {},
    thermalized_ {thermalized(
      graph,
      centerAtom,
      shape_,
      ranking_,
      Options::temperatureRegime
    )}
{}

/* Modification */
void AtomStereopermutator::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment && assignment.value() >= feasible_.indices.size()) {
    throw std::out_of_range("Supplied assignment index is out of range");
  }

  // Store new assignment
  assignmentOption_ = std::move(assignment);

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndex -> unsigned).
   */
  if(assignmentOption_) {
    shapePositionMap_ = siteToShapeVertexMap(
      abstract_.permutations.list.at(
        feasible_.indices.at(
          assignmentOption_.value()
        )
      ),
      abstract_.canonicalSites,
      ranking_.links
    );
  } else { // Wipe the map
    shapePositionMap_.clear();
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
          feasible_.indices,
          [&](const unsigned permutationIndex) -> unsigned {
            return abstract_.permutations.weights.at(permutationIndex);
          }
        ),
        engine
      )
    );
  }
}

void AtomStereopermutator::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // RankingInformation changes (lots of atom indices)
  ranking_.applyPermutation(permutation);

  // centerAtom_ must change
  centerAtom_ = permutation.at(centerAtom_);

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
#ifndef NDEBUG
  /* Check preconditions! */
  for(const auto& siteAtomList : boost::range::join(ranking_.sites, newRanking.sites)) {
    assert(
      std::is_sorted(
        std::begin(siteAtomList),
        std::end(siteAtomList)
      )
    );
  }
#endif

  /* There are a lot of changes that can occur prior to a call to propagate that
   * we need to correct. Multiple of these can apply at any given time!
   *
   * - Substituent count change
   *   - A new substituent may be added to a site or removed
   *   - Can happen at the same time as a site removal (three sites to two, one
   *     newly haptic)
   * - Ranking order change
   *   - Order of ranked sites may change
   *   - Implies a stereopermutation change!
   * - Site count change
   *   - Determined by comparing rankings
   *   - Implies a shape change
   *   - Need to make a map between old and new site indices (this can be
   *     complicated by substituent-level changes)
   * - Shape change
   *   - Triggered by a site count change only
   *   - Need to determine a map between old and new shape vertices
   *
   * Sketch of the following procedure
   * - Determine mapping between sites, taking care to account for site and
   *   substituent count changes
   * - Determine mapping of shape vertices between shapes if a shape change
   *   happens
   * - Use shape map of current assignment to place old site indices at old shape vertices
   * - Use mapping between shape vertices to enter the new shape
   * - Use mapping between site indices to place new sites at the new shape
   * - Make a stereopermutation from the new shape map and look for it in the
   *   new feasibles
   */

  // If nothing changes, nothing changes, and nothing is propagated
  if(newRanking == ranking_) {
    return boost::none;
  }

  const int siteCountChange = static_cast<int>(newRanking.sites.size() - ranking_.sites.size());
  // Complain if more than one site is changed either way
  if(std::abs(siteCountChange) > 1) {
    throw std::logic_error("propagateGraphChange should only ever handle a single site addition or removal at once");
  }

  auto siteMapping = SiteMapping::from(ranking_, newRanking);

  /* Decide the new shape */
  shapes::Shape newShape = shapeOption.value_or_eval(
    [&]() {
      if(siteCountChange == +1) {
        return up(shape_);
      }

      if(siteCountChange == 0) {
        return shape_;
      }

      assert(siteCountChange == -1);
      assert(siteMapping.changes);

      /* We can only figure out a mapping if the stereocenter is assigned
       * since otherwise cache_.symmetryPositionMap is empty and the shape
       * vertex being removed is a necessary argument to down.
       */
      if(assignmentOption_) {
        return down(
          shape_,
          shapePositionMap_.at(siteMapping.changes->site)
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
    centerAtom_,
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
   * ABCDEF. There will be multiple ways to split A to F!
   */
  if(assignmentOption_ && numStereopermutations() > 1) {
    const auto removedVertexOptional = temple::optionals::flatMap(
      siteMapping.changes,
      [this, siteCountChange](const SiteMapping::Changes changes) -> boost::optional<shapes::Vertex> {
        if(siteCountChange == -1) {
          return shapePositionMap_.at(changes.site);
        }

        return boost::none;
      }
    );

    const auto shapeMapping = temple::optionals::flatMap(
      shapes::getMapping(shape_, newShape, removedVertexOptional),
      std::bind(selectTransitionMapping, _1, Options::chiralStatePreservation)
    );

    const auto appliedShapeMapping = temple::optionals::map(
      shapeMapping,
      [&](const std::vector<shapes::Vertex>& vertexMap) {
        if(siteCountChange == +1) {
          return shapes::properties::applyIndexMapping(newShape, vertexMap);
        }

        return vertexMap;
      }
    );

    if(appliedShapeMapping) {
      const auto& shapeMappingValue = appliedShapeMapping.value();

      // Use shape map of current assignment to place old site indices at old shape vertices
      using InvertedMap = temple::StrongIndexFlatMap<shapes::Vertex, SiteIndex>;
      const InvertedMap oldSiteIndicesAtOldVertices = shapePositionMap_.invert();

      // Use shape mapping if available to transfer old vertices to new vertices
      const InvertedMap oldSiteIndicesAtNewVertices = InvertedMap(
        temple::map(
          temple::iota<shapes::Vertex>(shapes::size(newShape)),
          [&](const shapes::Vertex newVertex) -> SiteIndex {
            const shapes::Vertex oldVertex = shapeMappingValue.at(newVertex);

            if(siteCountChange == +1 && oldVertex == shapes::size(newShape) - 1)  {
              return siteMapping.changes->site;
            }

            return oldSiteIndicesAtOldVertices.at(oldVertex);
          }
        )
      );

      // Use site mapping to transfer old site indices to new site indices
      const InvertedMap sitesAtNewShapeVertices = InvertedMap(
        temple::map(
          oldSiteIndicesAtNewVertices,
          [&](const SiteIndex oldSite) -> SiteIndex {
            auto findIter = siteMapping.map.find(oldSite);

            if(findIter != std::end(siteMapping.map)) {
              return findIter->second;
            }

            assert(siteCountChange == +1);
            return oldSite;
          }
        )
      );

      assert(
        temple::all_of(
          sitesAtNewShapeVertices,
          [newShape](const SiteIndex i) { return i < shapes::size(newShape); }
        )
      );

      // Replace the site indices by their new ranking characters
      const auto newStereopermutationCharacters = stereopermutators::Abstract::makeStereopermutationCharacters(
        newAbstract.canonicalSites,
        newAbstract.symbolicCharacters,
        sitesAtNewShapeVertices
      );

      // Create a new assignment with those characters
      const auto trialStereopermutation = stereopermutation::Stereopermutation(
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
    std::move(ranking_),
    std::move(abstract_),
    std::move(feasible_),
    std::move(assignmentOption_)
  );

  // Overwrite the class state
  shape_ = newShape;
  ranking_ = std::move(newRanking);
  abstract_ = std::move(newAbstract);
  feasible_ = std::move(newFeasible);
  thermalized_ = thermalized(
    graph,
    centerAtom_,
    shape_,
    ranking_,
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
  assert(centerAtom_ != removedIndex);

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
  for(auto& equalPrioritySet : ranking_.substituentRanking) {
    for(auto& index : equalPrioritySet) {
      updateIndexInplace(index);
    }
  }

  for(auto& siteAtomList : ranking_.sites) {
    for(auto& atomIndex : siteAtomList) {
      updateIndexInplace(atomIndex);
    }
  }

  for(auto& link : ranking_.links) {
    link.cycleSequence = temple::map(link.cycleSequence, updateIndex);
  }
}

const stereopermutators::Abstract& AtomStereopermutator::Impl::getAbstract() const {
  return abstract_;
}

const stereopermutators::Feasible& AtomStereopermutator::Impl::getFeasible() const {
  return feasible_;
}

const RankingInformation& AtomStereopermutator::Impl::getRanking() const {
  return ranking_;
}

shapes::Shape AtomStereopermutator::Impl::getShape() const {
  return shape_;
}

const AtomStereopermutator::ShapeMap& AtomStereopermutator::Impl::getShapePositionMap() const {
  if(assignmentOption_ == boost::none) {
    throw std::logic_error(
      "The AtomStereopermutator is unassigned, sites are not assigned to "
      "symmetry positions"
    );
  }

  return shapePositionMap_;
}

void AtomStereopermutator::Impl::fit(
  const Graph& graph,
  const AngstromPositions& angstromWrapper
) {
  const unsigned S = shapes::size(shape_);
  assert(S == ranking_.sites.size());

  // Save stereopermutator state to return to if no fit is viable
  const shapes::Shape priorShape = shape_;
  const boost::optional<unsigned> priorStereopermutation  = assignmentOption_;

  // For all atoms making up a site, decide on the spatial average position
  Eigen::Matrix<double, 3, Eigen::Dynamic> sitePositions(3, S + 1);
  for(unsigned i = 0; i < S; ++i) {
    sitePositions.col(i) = cartesian::averagePosition(angstromWrapper.positions, ranking_.sites.at(i));
  }
  // Add the putative center
  sitePositions.col(S) = angstromWrapper.positions.row(centerAtom_);

  // Classify the shape and set it
  shapes::Shape fittedShape;
  std::vector<shapes::Vertex> matchingMapping;
  std::tie(fittedShape, matchingMapping) = classifyShape(sitePositions);
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
    SiteToShapeVertexMap {std::move(matchingMapping)},
    ranking_.links,
    abstract_.canonicalSites
  );

  auto soughtRotations = stereopermutation::generateAllRotations(soughtStereopermutation, fittedShape);

  /* Although the stereopermutations in abstract_ are sorted by their index of
   * permutation and we could potentially binary search all of our sought
   * rotations within that set, this linear search will never be time-critical
   * as long as we use continuous shape measure-based shape classification.
   */
  boost::optional<unsigned> foundStereopermutation;
  const unsigned A = feasible_.indices.size();
  for(unsigned a = 0; a < A; ++a) {
    const auto& feasiblePermutation = abstract_.permutations.list.at(
      feasible_.indices.at(a)
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
  const SiteIndex i,
  const SiteIndex j
) const {
  if(!assignmentOption_) {
    throw std::runtime_error("Stereopermutator is unassigned, angles are unknown!");
  }

  const unsigned S = shapes::size(shape_);
  if(i >= S || j >= S) {
    throw std::out_of_range("Site index is out of range");
  }

  return shapes::angleFunction(shape_)(
    shapePositionMap_.at(i),
    shapePositionMap_.at(j)
  );
}

boost::optional<unsigned> AtomStereopermutator::Impl::assigned() const {
  if(thermalized_) {
    return temple::optionals::map(assignmentOption_, [](unsigned /* a */) { return 0u; });
  }

  return assignmentOption_;
}

AtomIndex AtomStereopermutator::Impl::placement() const {
  return centerAtom_;
}

boost::optional<unsigned> AtomStereopermutator::Impl::indexOfPermutation() const {
  if(thermalized_) {
    return temple::optionals::map(assignmentOption_, [](unsigned /* a */) { return 0u; });
  }

  return temple::optionals::map(
    assignmentOption_,
    temple::functor::at(feasible_.indices)
  );
}

std::vector<AtomStereopermutator::MinimalChiralConstraint>
AtomStereopermutator::Impl::minimalChiralConstraints(const bool enforce) const {
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
    assignmentOption_
    && (numAssignments() > 1 || enforce)
  ) {
    /* Invert neighborSymmetryPositionMap_, we need a mapping of
     *  (vertex in shape) -> site index
     */
    const auto invertedShapeMap = shapePositionMap_.invert();

    return temple::map(
      shapes::tetrahedra(shape_),
      [&](const auto& tetrahedron) -> MinimalChiralConstraint {
        return temple::map_stl(
          tetrahedron,
          [&](const auto& shapeVertexOptional) -> boost::optional<SiteIndex> {
            return temple::optionals::map(
              shapeVertexOptional,
              temple::functor::at(invertedShapeMap)
            );
          }
        );
      }
    );
  }

  return {};
}

std::string AtomStereopermutator::Impl::info() const {
  std::string returnString = "A on "s
    + std::to_string(centerAtom_) + " ("s + shapes::name(shape_) +", "s;

  const auto& characters = abstract_.symbolicCharacters;
  std::copy(
    characters.begin(),
    characters.end(),
    std::back_inserter(returnString)
  );

  for(const auto& link : abstract_.selfReferentialLinks) {
    returnString += ", "s + characters.at(link.first) + "-"s + characters.at(link.second);
  }

  returnString += "): "s;

  if(assignmentOption_) {
    returnString += std::to_string(assignmentOption_.value());
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
    "A-"s + std::to_string(static_cast<unsigned>(shape_))
    + "-"s + std::to_string(numStereopermutations())
    + "-"s + (
      indexOfPermutation()
      ? std::to_string(indexOfPermutation().value())
      : "u"s
    )
  );
}

unsigned AtomStereopermutator::Impl::numAssignments() const {
  if(thermalized_) {
    return 1;
  }

  return feasible_.indices.size();
}

unsigned AtomStereopermutator::Impl::numStereopermutations() const {
  if(thermalized_) {
    return 1;
  }

  return abstract_.permutations.list.size();
}

void AtomStereopermutator::Impl::setShape(
  const shapes::Shape shape,
  const Graph& graph
) {
  if(shape_ == shape) {
    // If the symmetry doesn't actually change, then nothing does
    return;
  }

  shape_ = shape;

  abstract_ = stereopermutators::Abstract {
    ranking_,
    shape_
  };

  feasible_ = stereopermutators::Feasible {
    abstract_,
    shape_,
    centerAtom_,
    ranking_,
    graph
  };

  thermalized_ = thermalized(
    graph,
    centerAtom_,
    shape_,
    ranking_,
    Options::temperatureRegime
  );

  // Dis-assign the stereopermutator
  assign(boost::none);
}

} // namespace molassembler
} // namespace Scine
