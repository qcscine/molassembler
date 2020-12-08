/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutators/AtomStereopermutatorImpl.h"

#include "boost/range/join.hpp"
#include "Molassembler/Shapes/Properties.h"
#include "Molassembler/Shapes/PropertyCaching.h"
#include "Molassembler/Shapes/ContinuousMeasures.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "Molassembler/Stereopermutation/Manipulation.h"
#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Adaptors/Iota.h"
#include "Molassembler/Temple/Adaptors/Transform.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Functor.h"
#include "Molassembler/Temple/Optionals.h"
#include "Molassembler/Temple/OrderedPair.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/TinySet.h"
#include "Molassembler/Temple/constexpr/Math.h"
#include "Molassembler/Temple/constexpr/Numeric.h"

#include "Molassembler/Cycles.h"
#include "Molassembler/Detail/BuildTypeSwitch.h"
#include "Molassembler/Detail/Cartesian.h"
#include "Molassembler/DistanceGeometry/SpatialModel.h"
#include "Molassembler/DistanceGeometry/ValueBounds.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/Log.h"
#include "Molassembler/Modeling/CommonTrig.h"
#include "Molassembler/Stereopermutators/ShapeVertexMaps.h"
#include "Molassembler/Stereopermutators/RankingMapping.h"

#include "Utils/Geometry/ElementInfo.h"

using namespace std::string_literals;
using namespace std::placeholders;

namespace Scine {
namespace Molassembler {

namespace {

Shapes::Shape pickTransition(
  const Shapes::Shape shape,
  const unsigned T,
  const boost::optional<Shapes::Vertex>& removedVertexOptional
) {
  boost::optional<Shapes::Properties::ShapeTransitionGroup> bestTransition;
  std::vector<Shapes::Shape> propositions;

  auto replaceOrAdd = [&](
    Shapes::Shape newName,
    const Shapes::Properties::ShapeTransitionGroup& transitionGroup
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
  for(const Shapes::Shape propositionalShape : Shapes::allShapes) {
    if(Shapes::size(propositionalShape) != T) {
      continue;
    }

    if(auto transitionOptional = Shapes::getMapping(shape, propositionalShape, removedVertexOptional)) {
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
    return Shapes::Properties::mostSymmetric(T);
  }

  return Shapes::Properties::mostSymmetric(std::move(propositions));
}

std::pair<Shapes::Shape, std::vector<Shapes::Vertex>> classifyShape(const Eigen::Matrix<double, 3, Eigen::Dynamic>& sitePositions) {
  const unsigned S = sitePositions.cols() - 1;
  auto normalized = Shapes::Continuous::normalize(sitePositions);

  std::vector<Shapes::Shape> viableShapes;
  for(const Shapes::Shape shape : Shapes::allShapes) {
    if(Shapes::size(shape) == S) {
      viableShapes.push_back(shape);
    }
  }
  const unsigned shapesCount = viableShapes.size();
  std::vector<Shapes::Continuous::ShapeResult> shapeMeasureResults (shapesCount);
  std::vector<boost::optional<double>> randomCloudProbabilities (shapesCount);

#pragma omp parallel for
  for(unsigned i = 0; i < shapesCount; ++i) {
    const Shapes::Shape candidateShape = viableShapes[i];
    shapeMeasureResults[i] = Shapes::Continuous::shapeCentroidLast(normalized, candidateShape);
    // Shape classification for size 2 is better based on continuous shape measures themselves
    if(Shapes::size(candidateShape) > 2) {
      randomCloudProbabilities[i] = Shapes::Continuous::probabilityRandomCloud(shapeMeasureResults[i].measure, candidateShape);
    }
  }

  // Ensure centroids are mapped against one another
  assert(
    Temple::all_of(
      shapeMeasureResults,
      [](const auto& shapeResult) -> bool {
        assert(!shapeResult.mapping.empty());
        return shapeResult.mapping.back() == shapeResult.mapping.size() - 1;
      }
    )
  );

  // Prefer probabilities for comparison
  if(Temple::all_of(randomCloudProbabilities)) {
    const auto minElementIter = std::min_element(
      std::begin(randomCloudProbabilities),
      std::end(randomCloudProbabilities),
      [](const boost::optional<double>& a, const boost::optional<double>& b) -> bool {
        // We know all optionals are Somes from the previous all_of call
        return a.value() < b.value();
      }
    );
    const unsigned minimalShapeIndex = minElementIter - std::begin(randomCloudProbabilities);
    const Shapes::Shape minimalShape = viableShapes.at(minimalShapeIndex);
    return std::make_pair(
      minimalShape,
      std::move(shapeMeasureResults.at(minimalShapeIndex).mapping)
    );
  }

  // Fall back to minimal shape measure
  const auto minElementIter = std::min_element(
    std::begin(shapeMeasureResults),
    std::end(shapeMeasureResults),
    [](const auto& a, const auto& b) -> bool {
      return a.measure < b.measure;
    }
  );
  const unsigned minimalShapeIndex = minElementIter - std::begin(shapeMeasureResults);
  const Shapes::Shape minimalShape = viableShapes.at(minimalShapeIndex);
  return std::make_pair(
    minimalShape,
    std::move(minElementIter->mapping)
  );
}

} // namespace

Shapes::Shape AtomStereopermutator::Impl::up(const Shapes::Shape shape) {
  return pickTransition(shape, Shapes::size(shape) + 1, boost::none);
}

Shapes::Shape AtomStereopermutator::Impl::down(const Shapes::Shape shape, const Shapes::Vertex removedVertex) {
  return pickTransition(shape, Shapes::size(shape) - 1, removedVertex);
}

boost::optional<std::vector<Shapes::Vertex>> AtomStereopermutator::Impl::selectTransitionMapping(
  const Shapes::Properties::ShapeTransitionGroup& mappingsGroup,
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
    return Temple::Random::pick(mappingsGroup.indexMappings, randomnessEngine());
  }

  return boost::none;
}

bool AtomStereopermutator::Impl::thermalized(
  const Graph& graph,
  const AtomIndex centerAtom,
  const Shapes::Shape shape,
  const RankingInformation& ranking
) {
  if(Options::Thermalization::pyramidalInversion) {
    /* Nitrogen atom inversion */
    constexpr unsigned nitrogenZ = Utils::ElementInfo::Z(Utils::ElementType::N);
    const bool isNitrogenIsotope = Utils::ElementInfo::Z(graph.elementType(centerAtom)) == nitrogenZ;

    if(isNitrogenIsotope && shape == Shapes::Shape::VacantTetrahedron) {
      const bool inSmallCycle = Temple::any_of(
        ranking.links,
        [](const RankingInformation::Link& link) {
          return link.cycleSequence.size() <= 4;
        }
      );
      return !inSmallCycle;
    }
  }

  if(ranking.links.empty()) {
    if(Options::Thermalization::berryPseudorotation && shape == Shapes::Shape::TrigonalBipyramid) {
      return true;
    }

    if(Options::Thermalization::bartellMechanism && shape == Shapes::Shape::PentagonalBipyramid) {
      return true;
    }
  }

  return false;
}

/* Constructors */
AtomStereopermutator::Impl::Impl(
  const Graph& graph,
  // The symmetry of this Stereopermutator
  const Shapes::Shape shape,
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
      ranking_
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

void AtomStereopermutator::Impl::assign(std::vector<Shapes::Vertex> vertexMapping) {
  if(vertexMapping.size() != Shapes::size(shape_) + 1) {
    throw std::invalid_argument("Vertex mapping for assignment has wrong length");
  }

  /* Drop the centroid, making the remaining sequence viable as an index mapping
   * within the set of shape rotations, i.e. for use in stereopermutation
   * functionality
   */
  vertexMapping.pop_back();

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
    SiteToShapeVertexMap {vertexMapping},
    ranking_.links,
    abstract_.canonicalSites
  );

  auto soughtRotations = Stereopermutations::generateAllRotations(soughtStereopermutation, shape_);

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

  if(foundStereopermutation) {
    assignmentOption_ = *foundStereopermutation;
    shapePositionMap_ = ShapeMap {vertexMapping};
  }
}

void AtomStereopermutator::Impl::assignRandom(Random::Engine& engine) {
  const unsigned A = numAssignments();
  if(A == 0) {
    throw std::logic_error("Cannot randomly assign a stereopermutator without feasible stereopermutations");
  }

  if(A == 1) {
    assign(0);
  } else {
    assign(
      Temple::Random::pickDiscrete(
        // Map the feasible permutations onto their weights
        Temple::map(
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
  boost::optional<Shapes::Shape> shapeOption
) {
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
  const Shapes::Shape newShape = shapeOption.value_or_eval(
    [&]() {
      if(siteCountChange == +1) {
        return up(shape_);
      }

      if(siteCountChange == 0) {
        return shape_;
      }

      assert(siteCountChange == -1);
      assert(siteMapping.changedSite);

      /* We can only figure out a mapping if the stereocenter is assigned
       * since otherwise cache_.symmetryPositionMap is empty and the shape
       * vertex being removed is a necessary argument to down.
       */
      if(assignmentOption_) {
        return down(
          shape_,
          shapePositionMap_.at(siteMapping.changedSite.value())
        );
      }

      // Just return the most symmetric symmetry of the target size
      return Shapes::Properties::mostSymmetric(
        newRanking.sites.size()
      );
    }
  );

  /* Generate new assignments */
  Stereopermutators::Abstract newAbstract {
    newRanking,
    newShape
  };

  Stereopermutators::Feasible newFeasible {
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
    const auto removedVertexOptional = Temple::Optionals::flatMap(
      siteMapping.changedSite,
      [this, siteCountChange](const SiteIndex changedSite) -> boost::optional<Shapes::Vertex> {
        if(siteCountChange == -1) {
          return shapePositionMap_.at(changedSite);
        }

        return boost::none;
      }
    );

    const auto shapeMapping = Temple::Optionals::flatMap(
      Shapes::getMapping(shape_, newShape, removedVertexOptional),
      [](const auto& mapping) {
        return selectTransitionMapping(mapping, Options::chiralStatePreservation);
      }
    );

    const auto appliedShapeMapping = Temple::Optionals::map(
      shapeMapping,
      [&](const std::vector<Shapes::Vertex>& vertexMap) {
        if(siteCountChange == +1) {
          return Shapes::Properties::applyIndexMapping(newShape, vertexMap);
        }

        return vertexMap;
      }
    );

    if(appliedShapeMapping) {
      const auto& shapeMappingValue = appliedShapeMapping.value();

      // Use shape map of current assignment to place old site indices at old shape vertices
      using InvertedMap = Temple::StrongIndexFlatMap<Shapes::Vertex, SiteIndex>;
      const InvertedMap oldSiteIndicesAtOldVertices = shapePositionMap_.invert();

      // Use shape mapping if available to transfer old vertices to new vertices
      const InvertedMap oldSiteIndicesAtNewVertices = InvertedMap(
        Temple::map(
          Temple::iota<Shapes::Vertex>(Shapes::size(newShape)),
          [&](const Shapes::Vertex newVertex) -> SiteIndex {
            const Shapes::Vertex oldVertex = shapeMappingValue.at(newVertex);

            if(siteCountChange == +1 && oldVertex == Shapes::size(newShape) - 1)  {
              return siteMapping.changedSite.value();
            }

            return oldSiteIndicesAtOldVertices.at(oldVertex);
          }
        )
      );

      // Use site mapping to transfer old site indices to new site indices
      const InvertedMap sitesAtNewShapeVertices = InvertedMap(
        Temple::map(
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
        Temple::all_of(
          sitesAtNewShapeVertices,
          [newShape](const SiteIndex i) { return i < Shapes::size(newShape); }
        )
      );

      // Replace the site indices by their new ranking characters
      const auto newStereopermutationCharacters = Stereopermutators::Abstract::makeStereopermutationCharacters(
        newAbstract.canonicalSites,
        newAbstract.symbolicCharacters,
        sitesAtNewShapeVertices
      );

      // Create a new assignment with those characters
      const auto trialStereopermutation = Stereopermutations::Stereopermutation(
        newStereopermutationCharacters,
        newAbstract.selfReferentialLinks
      );

      // Generate all rotations of this trial assignment
      auto allTrialRotations = Stereopermutations::generateAllRotations(trialStereopermutation, newShape);

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
  auto newAssignmentOption = Temple::Optionals::flatMap(
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

  if(!newAssignmentOption && newFeasible.indices.size() == 1) {
    newAssignmentOption = 0;
  }

  // Extract old state from the class
  auto oldStateTuple = std::make_tuple(
    std::move(ranking_),
    std::move(abstract_),
    std::move(feasible_),
    std::move(shapePositionMap_)
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
    ranking_
  );
  assign(newAssignmentOption);

  return {std::move(oldStateTuple)};
}

void AtomStereopermutator::Impl::propagateVertexRemoval(const AtomIndex removedIndex) {
  /* This function replaces any occurrences of the atom index that is being
   * removed in the global state with a placeholder of the same type and updates
   * any invalidated atom indices.
   */

  /* If the central atom is being removed, drop this stereopermutator
   * beforehand in caller. This would just be unnecessary work.
   */
  assert(centerAtom_ != removedIndex);

  auto updateIndex = [&removedIndex](const AtomIndex index) -> AtomIndex {
    if(index > removedIndex) {
      return index - 1;
    }

    if(index == removedIndex) {
      return PrivateGraph::removalPlaceholder;
    }

    return index;
  };

  centerAtom_ = updateIndex(centerAtom_);

  /* Update indices in RankingInformation */
  for(auto& equalPrioritySet : ranking_.substituentRanking) {
    for(auto& index : equalPrioritySet) {
      index = updateIndex(index);
    }
  }

  for(auto& siteAtomList : ranking_.sites) {
    for(auto& atomIndex : siteAtomList) {
      atomIndex = updateIndex(atomIndex);
    }
  }

  for(auto& link : ranking_.links) {
    link.cycleSequence = Temple::map(link.cycleSequence, updateIndex);
  }
}

const Stereopermutators::Abstract& AtomStereopermutator::Impl::getAbstract() const {
  return abstract_;
}

const Stereopermutators::Feasible& AtomStereopermutator::Impl::getFeasible() const {
  return feasible_;
}

const RankingInformation& AtomStereopermutator::Impl::getRanking() const {
  return ranking_;
}

Shapes::Shape AtomStereopermutator::Impl::getShape() const {
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

boost::optional<AtomStereopermutator::ShapeMap>
AtomStereopermutator::Impl::fit(
  const Graph& graph,
  const AngstromPositions& angstromWrapper
) {
  const unsigned S = Shapes::size(shape_);
  assert(S == ranking_.sites.size());

  // Save stereopermutator state to return to if no fit is viable
  const Shapes::Shape priorShape = shape_;
  const boost::optional<unsigned> priorStereopermutation  = assignmentOption_;

  // For all atoms making up a site, decide on the spatial average position
  Eigen::Matrix<double, 3, Eigen::Dynamic> sitePositions(3, S + 1);
  for(unsigned i = 0; i < S; ++i) {
    sitePositions.col(i) = Cartesian::averagePosition(angstromWrapper.positions, ranking_.sites.at(i));
  }
  // Add the putative center
  sitePositions.col(S) = angstromWrapper.positions.row(centerAtom_);

  // Classify the shape and set it
  Shapes::Shape fittedShape;
  std::vector<Shapes::Vertex> matchingMapping;
  std::tie(fittedShape, matchingMapping) = classifyShape(sitePositions);
  setShape(fittedShape, graph);
  assign(matchingMapping);

  // Return to prior state if no feasible stereopermutations found
  if(assignmentOption_ == boost::none) {
    setShape(priorShape, graph);
    assign(priorStereopermutation);
    return boost::none;
  }

  matchingMapping.pop_back();
  return ShapeMap(matchingMapping);
}

/* Information */
double AtomStereopermutator::Impl::angle(
  const SiteIndex i,
  const SiteIndex j
) const {
  if(!assignmentOption_) {
    throw std::runtime_error("Stereopermutator is unassigned, angles are unknown!");
  }

  const unsigned S = Shapes::size(shape_);
  if(i >= S || j >= S) {
    throw std::out_of_range("Site index is out of range");
  }

  return Shapes::angleFunction(shape_)(
    shapePositionMap_.at(i),
    shapePositionMap_.at(j)
  );
}

boost::optional<unsigned> AtomStereopermutator::Impl::assigned() const {
  if(thermalized_) {
    return Temple::Optionals::map(assignmentOption_, [](unsigned /* a */) { return 0u; });
  }

  return assignmentOption_;
}

AtomIndex AtomStereopermutator::Impl::placement() const {
  return centerAtom_;
}

boost::optional<unsigned> AtomStereopermutator::Impl::indexOfPermutation() const {
  if(thermalized_) {
    return Temple::Optionals::map(assignmentOption_, [](unsigned /* a */) { return 0u; });
  }

  return Temple::Optionals::map(
    assignmentOption_,
    Temple::Functor::at(feasible_.indices)
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

    return Temple::map(
      Shapes::tetrahedra(shape_),
      [&](const auto& tetrahedron) -> MinimalChiralConstraint {
        return Temple::map(
          tetrahedron,
          [&](const auto& shapeVertexOptional) -> boost::optional<SiteIndex> {
            return Temple::Optionals::map(
              shapeVertexOptional,
              Temple::Functor::at(invertedShapeMap)
            );
          }
        );
      }
    );
  }

  return {};
}

std::string AtomStereopermutator::Impl::info() const {
  std::string returnString = std::to_string(centerAtom_) + ": "s + Shapes::name(shape_) +", "s;

  const auto& characters = abstract_.symbolicCharacters;
  std::copy(
    characters.begin(),
    characters.end(),
    std::back_inserter(returnString)
  );

  for(const auto& link : abstract_.selfReferentialLinks) {
    returnString += ", "s + characters.at(link.first) + "-"s + characters.at(link.second);
  }

  returnString += ", "s;

  if(assignmentOption_) {
    returnString += std::to_string(assignmentOption_.value());
  } else {
    returnString += "u";
  }

  const unsigned A = numAssignments();
  returnString += " ("s + std::to_string(A);

  const unsigned P = numStereopermutations();
  if(P != A) {
    returnString += ", "s + std::to_string(P);
  }

  returnString += ")";

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

std::vector<std::vector<SiteIndex>> AtomStereopermutator::Impl::siteGroups() const {
  if(!assignmentOption_) {
    throw std::logic_error("Stereopermutator is not assigned");
  }

  return Temple::map(
    Shapes::Properties::positionGroups(shape_),
    [this](const std::vector<Shapes::Vertex>& interconvertibleVertices) {
      return Temple::map(
        interconvertibleVertices,
        [this](const Shapes::Vertex v) -> SiteIndex {
          return shapePositionMap_.indexOf(v);
        }
      );
    }
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
  const Shapes::Shape shape,
  const Graph& graph
) {
  if(shape_ == shape) {
    // If the symmetry doesn't actually change, then nothing does
    return;
  }

  shape_ = shape;

  abstract_ = Stereopermutators::Abstract {
    ranking_,
    shape_
  };

  feasible_ = Stereopermutators::Feasible {
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
    ranking_
  );

  // Dis-assign the stereopermutator
  assign(boost::none);
}

} // namespace Molassembler
} // namespace Scine
