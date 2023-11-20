/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Molecule/MoleculeImpl.h"

#include "Molassembler/Shapes/constexpr/Properties.h"
#include "boost/functional/hash.hpp"
#include "boost/graph/graph_utility.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/isomorphism.hpp"

#include "Molassembler/Cycles.h"
#include "Molassembler/Detail/Cartesian.h"
#include "Molassembler/Graph/Canonicalization.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Modeling/CommonTrig.h"
#include "Molassembler/Modeling/ShapeInference.h"
#include "Molassembler/Molecule/AtomEnvironmentHash.h"
#include "Molassembler/Molecule/MolGraphWriter.h"
#include "Molassembler/Molecule/RankingTree.h"
#include "Molassembler/Options.h"
#include "Molassembler/Stereopermutation/Composites.h"
#include "Molassembler/Stereopermutators/AbstractPermutations.h"
#include "Molassembler/Stereopermutators/FeasiblePermutations.h"

#include "Molassembler/Temple/Optionals.h"

namespace Scine {
namespace Molassembler {

Utils::AtomCollection Molecule::Impl::applyCanonicalizationMap(
  const std::vector<AtomIndex>& canonicalizationIndexMap,
  const Utils::AtomCollection& atomCollection) {
  const int N = atomCollection.size();
  Utils::AtomCollection permuted(N);
  for(int i = 0; i < N; ++i) {
    int newIndex = canonicalizationIndexMap.at(i);
    permuted.setElement(newIndex, atomCollection.getElement(i));
    permuted.setPosition(newIndex, atomCollection.getPosition(i));
  }

  return permuted;
}

boost::optional<AtomStereopermutator> Molecule::Impl::makePermutator(
  AtomIndex candidateIndex,
  const StereopermutatorList& stereopermutators,
  const boost::optional<AngstromPositions>& maybePositions,
  const boost::optional<SubstitutionsGenerator::SubstitutionMap>& maybeSubstitutions) const {
  // If there is already an atom stereopermutator on this index, stop
  if(stereopermutators.option(candidateIndex)) {
    return boost::none;
  }

  RankingInformation ranking = rankPriority(candidateIndex, {}, maybePositions);

  // Only non-terminal atoms may have permutators
  const unsigned S = ranking.sites.size();
  if(S <= 1) {
    return boost::none;
  }

  const Shapes::Shape shape = [&]() {
    if(maybePositions) {
      /* If positions are present, then the shape will be fitted, not inferred,
       * so the shape in construction is of little consequence
       */
      return ShapeInference::firstOfSize(S);
    }

    return inferShape(candidateIndex, ranking).value_or_eval([&]() {
      return ShapeInference::firstOfSize(S);
    });
  }();

  const Stereopermutators::Feasible::Functor feasibilityFunctor{adjacencies_};
  const auto thermalizationFunctor
    = AtomStereopermutator::thermalizationFunctor(adjacencies_);

  AtomStereopermutator permutator{
    candidateIndex, shape, std::move(ranking), feasibilityFunctor, thermalizationFunctor};

  if(maybePositions) {
    const auto substitutions
      = Temple::Optionals::flatMap(
          maybeSubstitutions,
          [&](const auto& map) {
            return Temple::Optionals::mapFind(map, candidateIndex);
          })
          .value_or(SubstitutionsGenerator::SubstitutionList{});

    const auto siteCentroids
      = permutator.sitePositions(maybePositions.value(), substitutions);
    permutator.fit(siteCentroids, feasibilityFunctor, thermalizationFunctor);
  } else if(permutator.numAssignments() == 1) {
    permutator.assign(0);
  }
  return permutator;
}

boost::optional<BondStereopermutator> Molecule::Impl::makePermutator(
  const BondIndex& bond,
  const StereopermutatorList& stereopermutators,
  const boost::optional<AngstromPositions>& maybePositions,
  const boost::optional<SubstitutionsGenerator::SubstitutionMap>& maybeSubstitutions,
  const BondStereopermutator::Alignment alignment) const {
  // If there is already a bond stereopermutator on this edge, stop
  if(stereopermutators.option(bond)) {
    return boost::none;
  }

  const auto permutatorOptions = Temple::map(
    bond, [&](const AtomIndex side) { return stereopermutators.option(side); });

  const auto existAndAssigned
    = Temple::map(permutatorOptions, [](const auto& permutatorOption) -> bool {
        return permutatorOption && permutatorOption->assigned() != boost::none;
      });

  if(!existAndAssigned.first || !existAndAssigned.second) {
    return boost::none;
  }

  // If one of the permutators is considered thermalized we must not construct
  // a bond stereopermutator because at least one atom is considered to occupy
  // multiple permutations at once.
  if(permutatorOptions.first->thermalized() || permutatorOptions.second->thermalized()) {
    return boost::none;
  }

  // Construct a Stereopermutator here
  BondStereopermutator permutator{adjacencies_.inner(), stereopermutators, bond, alignment};

  // Default-assign single-assignment bond stereopermutators
  if(maybePositions) {
    const auto fittingRefs = Temple::map(
      permutatorOptions,
      [&](const auto& permutatorOption) -> BondStereopermutator::FittingReferences {
        return {permutatorOption.value(), permutatorOption->getShapePositionMap()};
      });

    const auto sitePositions = permutator.composite().orientations().map(
      [&](const auto& orientation) -> AtomStereopermutator::SiteCentroids {
        const auto substitutions
          = Temple::Optionals::flatMap(
              maybeSubstitutions,
              [&](const auto& map) {
                return Temple::Optionals::mapFind(map, orientation.identifier);
              })
              .value_or(SubstitutionsGenerator::SubstitutionList{});

        return stereopermutators.at(orientation.identifier)
          .sitePositions(maybePositions.value(), substitutions);
      });

    permutator.fit(sitePositions, fittingRefs);
  } else if(permutator.numAssignments() == 1) {
    permutator.assign(0);
  }

  return permutator;
}

StereopermutatorList Molecule::Impl::detectStereopermutators_() const {
  StereopermutatorList stereopermutatorList;

#ifdef _OPENMP
  /* Ensure inner's properties are populated to avoid data races in its mutable
   * members
   */
  adjacencies_.inner().populateProperties();
#endif

  for(const AtomIndex candidate : graph().atoms()) {
    if(auto permutator = makePermutator(candidate, stereopermutatorList)) {
      stereopermutatorList.add(std::move(permutator.value()));
    }
  }

  for(BondIndex bond : graph().bonds()) {
    if(!isGraphBasedBondStereopermutatorCandidate_(graph().bondType(bond))) {
      continue;
    }

    if(auto permutator = makePermutator(bond, stereopermutatorList)) {
      if(permutator->numStereopermutations() > 1) {
        stereopermutatorList.add(std::move(permutator.value()));
      }
    }
  }

  return stereopermutatorList;
}

void Molecule::Impl::ensureModelInvariants_() const {
  if(graph().inner().connectedComponents() > 1) {
    throw std::logic_error("Molecules must be a single connected component. "
                           "The supplied graph has multiple");
  }

  if(graph().V() < 1) {
    throw std::logic_error("Molecules must consist of at least one atom!");
  }
}

bool Molecule::Impl::isValidIndex_(const AtomIndex index) const {
  return index < graph().V();
}

bool Molecule::Impl::isGraphBasedBondStereopermutatorCandidate_(BondType bondType) {
  return (
    bondType == BondType::Double || bondType == BondType::Triple
    || bondType == BondType::Quadruple || bondType == BondType::Quintuple
    || bondType == BondType::Sextuple);
}

void Molecule::Impl::propagateGraphChange_() {
  /* Two cases: If the StereopermutatorList is empty, we can just use detect to
   * find any new stereopermutators in the Molecule.
   */
  if(stereopermutators_.empty()) {
    stereopermutators_ = detectStereopermutators_();
    return;
  }

  /* In the other case, we have to recheck absolutely everywhere because of the
   * IUPAC sequence rules. If a stereopermutator's ranking was affected and the
   * stereopermutator has a set assignment, we then need to find the assignment
   * that the previous ranking represented spatially in the new set of
   * assignments and assign the stereopermutator to that.
   */

  GraphAlgorithms::updateEtaBonds(adjacencies_.inner());

  const PrivateGraph& inner = adjacencies_.inner();

  for(const PrivateGraph::Vertex vertex : inner.vertices()) {
    auto stereopermutatorOption = stereopermutators_.option(vertex);
    RankingInformation localRanking = rankPriority(vertex);

    if(stereopermutatorOption) {
      // The atom has become terminal
      if(localRanking.sites.size() <= 1) {
        stereopermutators_.remove(vertex);
        continue;
      }

      // Has the ranking changed?
      if(localRanking == stereopermutatorOption->getRanking()) {
        continue;
      }

      // Are there adjacent bond stereopermutators?
      std::vector<BondIndex> adjacentStereopermutatorBonds;
      for(BondIndex bond : adjacencies_.bonds(vertex)) {
        if(stereopermutators_.option(bond)) {
          adjacentStereopermutatorBonds.push_back(std::move(bond));
        }
      }

      // Suggest a shape if desired
      boost::optional<Shapes::Shape> newShapeOption;
      if(Options::shapeTransition == ShapeTransition::PrioritizeInferenceFromGraph) {
        newShapeOption = inferShape(vertex, localRanking);
      }

      // Propagate the state
      auto oldAtomStereopermutatorStateOption = stereopermutatorOption->propagate(
        std::move(localRanking),
        newShapeOption,
        Stereopermutators::Feasible::Functor(adjacencies_),
        AtomStereopermutator::thermalizationFunctor(adjacencies_));

      /* If the chiral state for this atom stereopermutator was not successfully
       * propagated or it is now unassigned, then bond stereopermutators sharing
       * this atom stereopermutator must be removed. Bond stereopermutators can
       * only be undetermined if its constituting atom stereopermutators are
       * assigned.
       */
      if(!stereopermutatorOption->assigned()) {
        for(const BondIndex& bond : adjacentStereopermutatorBonds) {
          stereopermutators_.remove(bond);
        }

        continue;
      }

      /* If the chiral state for this atom stereopermutator was successfully
       * propagated and/or the permutator could be default-assigned, we can also
       * propagate adjacent BondStereopermutators.
       *
       * TODO we may have to keep track if assignments change within the
       * propagated bondstereopermutators, or if any bond stereopermutators
       * are removed, since this may cause another re-rank!
       */
      if(oldAtomStereopermutatorStateOption) {
        for(const BondIndex& bond : adjacentStereopermutatorBonds) {
          // stereopermutators_.remove(bond);
          // stereopermutators_.add(BondStereopermutator(inner, stereopermutators_, bond));
          stereopermutators_.at(bond).propagateGraphChange(
            *oldAtomStereopermutatorStateOption,
            *stereopermutatorOption,
            inner,
            stereopermutators_
          );
        }
      }
    } else if(auto permutator = makePermutator(vertex, stereopermutators_)) {
      // There is no atom stereopermutator on this vertex, so try to add one
      stereopermutators_.add(std::move(*permutator));
    }
  }

  // Look for new bond stereopermutators
  for(const BondIndex bond : graph().bonds()) {
    if(
      stereopermutators_.option(bond)
      || !isGraphBasedBondStereopermutatorCandidate_(graph().bondType(bond))) {
      continue;
    }

    if(auto permutator = makePermutator(bond, stereopermutators_)) {
      if(permutator->numStereopermutations() > 1) {
        stereopermutators_.add(std::move(*permutator));
      }
    }
  }
}

/* Public members */
/* Constructors */
Molecule::Impl::Impl() noexcept :
  Impl(Utils::ElementType::H, Utils::ElementType::H, BondType::Single) {
}

Molecule::Impl::Impl(const Utils::ElementType element) noexcept {
  PrivateGraph& inner = adjacencies_.inner();
  inner.addVertex(element);
}

Molecule::Impl::Impl(const Utils::ElementType a, const Utils::ElementType b, const BondType bondType) noexcept {
  // update adjacencies_
  PrivateGraph& inner = adjacencies_.inner();
  PrivateGraph::Vertex i = inner.addVertex(a);
  PrivateGraph::Vertex j = inner.addVertex(b);
  inner.addEdge(i, j, bondType);
}

Molecule::Impl::Impl(Graph graph) : adjacencies_(std::move(graph)) {
  // Initialization
  GraphAlgorithms::updateEtaBonds(adjacencies_.inner());
  stereopermutators_ = detectStereopermutators_();
  ensureModelInvariants_();
}

Molecule::Impl::Impl(
  Graph graph,
  const AngstromPositions& positions,
  const boost::optional<std::vector<BondIndex>>& bondStereopermutatorCandidatesOptional,
  const boost::optional<PeriodicBoundaryDuplicates>& periodics) :
  adjacencies_(std::move(graph)) {
  const auto maybeSubstitutions
    = Temple::Optionals::map(periodics, [&](const PeriodicBoundaryDuplicates& some) {
        return SubstitutionsGenerator::removeGhosts(adjacencies_, some);
      });

  const auto maybeUninteresting = Temple::Optionals::map(
    periodics, [](const auto& p) -> const auto& { return p.uninterestingAtoms; });

  GraphAlgorithms::updateEtaBonds(adjacencies_.inner());
  stereopermutators_ = inferStereopermutatorsFromPositions(
    positions, bondStereopermutatorCandidatesOptional, maybeUninteresting, maybeSubstitutions);
  ensureModelInvariants_();
}

Molecule::Impl::Impl(
  Graph graph,
  StereopermutatorList stereopermutators,
  boost::optional<AtomEnvironmentComponents> canonicalComponentsOption) :
  adjacencies_(std::move(graph)),
  stereopermutators_(std::move(stereopermutators)),
  canonicalComponentsOption_(std::move(canonicalComponentsOption)) {
  // Initialization
  ensureModelInvariants_();
}

/* Modifiers */
AtomIndex Molecule::Impl::addAtom(
  const Utils::ElementType elementType,
  const AtomIndex adjacentTo,
  const BondType bondType) {
  if(!isValidIndex_(adjacentTo)) {
    throw std::out_of_range(
      "Molecule::addAtom: Supplied atom index is invalid!");
  }

  const AtomIndex index = adjacencies_.inner().addVertex(elementType);
  addBond(index, adjacentTo, bondType);
  /* addBond handles the stereopermutator update on adjacentTo and also
   * performs a full-molecule rerank propagation.
   */

  return index;
}

BondIndex Molecule::Impl::addBond(const AtomIndex a, const AtomIndex b, const BondType bondType) {
  if(!isValidIndex_(a) || !isValidIndex_(b)) {
    throw std::out_of_range("Molecule::addBond: A supplied index is invalid!");
  }

  if(a == b) {
    throw std::logic_error(
      "Molecule::addBond: Cannot add a bond between identical indices!");
  }

  PrivateGraph& inner = adjacencies_.inner();

  inner.addEdge(a, b, bondType);

  auto notifySubstituentAddition = [this](const AtomIndex toIndex) {
    /*! @todo Remove any BondStereopermutators on adjacent edges of toIndex (no
     * substituent addition/removal propagation possible yet)
     */
    for(const BondIndex& adjacentEdge : adjacencies_.bonds(toIndex)) {
      stereopermutators_.remove(adjacentEdge);
    }
  };

  notifySubstituentAddition(a);
  notifySubstituentAddition(b);

  propagateGraphChange_();
  canonicalComponentsOption_ = boost::none;

  return BondIndex{a, b};
}

const BondStereopermutator& Molecule::Impl::addPermutator(
  const BondIndex& bond,
  const BondStereopermutator::Alignment alignment) {
  auto permutator
    = makePermutator(bond, stereopermutators_, boost::none, boost::none, alignment);
  if(!permutator) {
    throw std::logic_error("Violated preconditions for permutator addition");
  }
  stereopermutators_.add(std::move(permutator.value()));
  propagateGraphChange_();
  return stereopermutators_.at(bond);
}

void Molecule::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
  adjacencies_.inner().applyPermutation(permutation);
  stereopermutators_.applyPermutation(permutation);
  canonicalComponentsOption_ = boost::none;
}

void Molecule::Impl::assignStereopermutator(
  const AtomIndex a,
  const boost::optional<unsigned>& assignmentOption) {
  if(!isValidIndex_(a)) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  auto stereopermutatorOption = stereopermutators_.option(a);

  if(!stereopermutatorOption) {
    throw std::out_of_range("No stereopermutator at this index!");
  }

  if(assignmentOption && assignmentOption.value() >= stereopermutatorOption->numAssignments()) {
    throw std::out_of_range("Invalid assignment index!");
  }

  // No need to rerank and remove canonical components if nothing changes
  if(assignmentOption != stereopermutatorOption->assigned()) {
    stereopermutatorOption->assign(assignmentOption);

    // A reassignment can change ranking! See the RankingTree tests
    propagateGraphChange_();
    canonicalComponentsOption_ = boost::none;
  }
}

void Molecule::Impl::assignStereopermutator(
  const BondIndex& edge,
  const boost::optional<unsigned>& assignmentOption) {
  if(!isValidIndex_(edge.first) || !isValidIndex_(edge.second)) {
    throw std::out_of_range("Supplied bond atom indices are invalid!");
  }

  auto stereopermutatorOption = stereopermutators_.option(edge);

  if(!stereopermutatorOption) {
    throw std::out_of_range(
      "assignStereopermutator: No stereopermutator at this bond!");
  }

  if(assignmentOption && assignmentOption.value() >= stereopermutatorOption->numAssignments()) {
    throw std::out_of_range(
      "assignStereopermutator: Invalid assignment index!");
  }

  // No need to rerank and remove canonical components if nothing changes
  if(assignmentOption != stereopermutatorOption->assigned()) {
    stereopermutatorOption->assign(assignmentOption);

    // A reassignment can change ranking! See the RankingTree tests
    propagateGraphChange_();
    canonicalComponentsOption_ = boost::none;
  }
}

void Molecule::Impl::assignStereopermutatorRandomly(const AtomIndex a, Random::Engine& engine) {
  if(!isValidIndex_(a)) {
    throw std::out_of_range(
      "Molecule::assignStereopermutatorRandomly: Supplied index is invalid!");
  }

  auto stereopermutatorOption = stereopermutators_.option(a);

  if(!stereopermutatorOption) {
    throw std::out_of_range(
      "assignStereopermutatorRandomly: No stereopermutator at this index!");
  }

  stereopermutatorOption->assignRandom(engine);

  // A reassignment can change ranking! See the RankingTree tests
  propagateGraphChange_();
  canonicalComponentsOption_ = boost::none;
}

void Molecule::Impl::assignStereopermutatorRandomly(const BondIndex& e, Random::Engine& engine) {
  auto stereopermutatorOption = stereopermutators_.option(e);

  if(!stereopermutatorOption) {
    throw std::out_of_range(
      "assignStereopermutatorRandomly: No stereopermutator at this edge!");
  }

  stereopermutatorOption->assignRandom(engine);

  // A reassignment can change ranking! See the RankingTree tests
  propagateGraphChange_();
  canonicalComponentsOption_ = boost::none;
}

std::vector<AtomIndex>
Molecule::Impl::canonicalize(const AtomEnvironmentComponents componentBitmask) {
  // Generate hashes according to the passed bitmask
  auto vertexHashes
    = Hashes::generate(graph().inner(), stereopermutators(), componentBitmask);

  // Get a canonical labelling
  auto labelMap = canonicalAutomorphism(graph().inner(), vertexHashes);

  struct Inverse {
    std::vector<AtomIndex> permutation;

    inline explicit Inverse(const std::vector<int>& ante) {
      unsigned size = ante.size();
      permutation.resize(size);
      for(AtomIndex i = 0; i < size; ++i) {
        permutation.at(ante.at(i)) = i;
      }
    }

    inline AtomIndex operator()(const AtomIndex i) const {
      return permutation.at(i);
    }
  };

  Inverse inverse(labelMap);

  // Apply the labelling change to the graph and all stereopermutators
  applyPermutation(inverse.permutation);

  // Set the underlying canonical components
  canonicalComponentsOption_ = componentBitmask;

  return inverse.permutation;
}

void Molecule::Impl::removeAtom(const AtomIndex a) {
  if(!isValidIndex_(a)) {
    throw std::out_of_range("Molecule::removeAtom: Supplied index is invalid!");
  }

  if(!graph().canRemove(a)) {
    throw std::logic_error("Removing this atom disconnects the graph!");
  }

  PrivateGraph& inner = adjacencies_.inner();

  std::vector<AtomIndex> previouslyAdjacentVertices;
  auto adjacencyIterators = inner.adjacents(a);
  std::copy(
    adjacencyIterators.first,
    adjacencyIterators.second,
    std::back_inserter(previouslyAdjacentVertices));

  // Any adjacent bond stereopermutators have to be dropped
  for(const BondIndex& adjacentEdge : adjacencies_.bonds(a)) {
    stereopermutators_.remove(adjacentEdge);
  }

  // Remove all edges to and from this vertex
  inner.clearVertex(a);

  // Update Eta bonds
  GraphAlgorithms::updateEtaBonds(adjacencies_.inner());

  // Any stereopermutator on this index must be dropped
  stereopermutators_.remove(a);

  // Remove the vertex itself
  inner.removeVertex(a);

  /* Removing the vertex invalidates some vertex descriptors, which are used
   * liberally in the stereopermutator classes' state. We have to correct
   * all vertex descriptors across all permutators to ensure that chiral state
   * is modeled correctly.
   */
  stereopermutators_.propagateVertexRemoval(a);

  /* call removeSubstituent on all adjacent stereopermutators, with
   * removalPlaceholder as the 'which' parameter, which is what
   * propagateVertexRemoval replaces the removed index with in the
   * stereopermutators' internal state
   */

  for(AtomIndex indexToUpdate : previouslyAdjacentVertices) {
    if(indexToUpdate > a) {
      --indexToUpdate;
    }

    auto stereopermutatorOption = stereopermutators_.option(indexToUpdate);
    if(!stereopermutatorOption) {
      /* If there's no atom stereopermutator on this vertex, there's nothing
       * to do
       */
      continue;
    }

    RankingInformation localRanking = rankPriority(indexToUpdate);

    /* If index becomes terminal, drop the stereopermutator and any adjoining
     * bond stereopermutators
     */

    if(localRanking.sites.size() <= 1) {
      stereopermutators_.remove(indexToUpdate);
      for(const BondIndex edge : adjacencies_.bonds(indexToUpdate)) {
        stereopermutators_.remove(edge);
      }
      continue;
    }

    // Generate a target shape for the atom stereopermutator
    boost::optional<Shapes::Shape> newShapeOption;
    // if(Options::shapeTransition ==
    // ShapeTransition::PrioritizeInferenceFromGraph) {
    newShapeOption = inferShape(indexToUpdate, localRanking);
    // }

    const auto propagatedState = stereopermutatorOption->propagate(
      std::move(localRanking),
      newShapeOption,
      Stereopermutators::Feasible::Functor(adjacencies_),
      AtomStereopermutator::thermalizationFunctor(adjacencies_));

    for(const BondIndex edge : adjacencies_.bonds(indexToUpdate)) {
      if(propagatedState && stereopermutators_.option(edge)) {
        if(
          !(stereopermutators_.option(edge.first)->assigned())
          || !(stereopermutators_.option(edge.second)->assigned())) {
          continue;
        }
        stereopermutators_.remove(edge);

        stereopermutators_.add(BondStereopermutator(inner, stereopermutators_, edge));

      } else if(stereopermutators_.option(edge)) {
        stereopermutators_.remove(edge);
      }
    }
  }

  propagateGraphChange_();

  canonicalComponentsOption_ = boost::none;
}

void Molecule::Impl::removeBond(const AtomIndex a, const AtomIndex b) {
  if(!isValidIndex_(a) || !isValidIndex_(b)) {
    throw std::out_of_range("Molecule::removeBond: Supplied index is invalid!");
  }

  PrivateGraph& inner = adjacencies_.inner();

  auto edgeOption = inner.edgeOption(a, b);

  if(!edgeOption) {
    throw std::out_of_range("That bond does not exist!");
  }

  PrivateGraph::Edge edgeToRemove = edgeOption.value();

  if(!inner.canRemove(edgeToRemove)) {
    throw std::logic_error(
      "Removing this bond separates the molecule into two pieces!");
  }

  /* If there is an BondStereopermutator on this edge, we have to drop it
   * explicitly, since propagateGraphChange_ cannot iterate over a now-removed
   * edge.
   */
  stereopermutators_.remove(BondIndex{a, b});

  // Remove the edge
  inner.removeEdge(edgeToRemove);

  // Notify all immediately adjacent stereopermutators of the removal
  auto notifyRemoval = [&](const AtomIndex indexToUpdate) {
    if(auto stereopermutatorOption = stereopermutators_.option(indexToUpdate)) {
      auto localRanking = rankPriority(indexToUpdate);

      // In case the central atom becomes terminal, just drop the stereopermutator
      if(localRanking.sites.size() <= 1) {
        stereopermutators_.remove(indexToUpdate);
        return;
      }

      boost::optional<Shapes::Shape> newShapeOption;
      if(Options::shapeTransition == ShapeTransition::PrioritizeInferenceFromGraph) {
        newShapeOption = inferShape(indexToUpdate, localRanking);
      }

      stereopermutatorOption->propagate(
        std::move(localRanking),
        newShapeOption,
        Stereopermutators::Feasible::Functor(adjacencies_),
        AtomStereopermutator::thermalizationFunctor(adjacencies_));
    }

    //! @todo propagation
    /*if(stereopermutators_.involving(indexToUpdate)) {
      if(stereopermutators_.at(indexToUpdate)->type() ==
    Stereopermutators::Type::AtomStereopermutator) { } else {
        std::dynamic_pointer_cast<Stereopermutators::BondStereopermutator>(
          stereopermutators_.at(indexToUpdate)
        )->removeSubstituent(
          indexToUpdate,
          removedIndex
        );
      }
    }*/
  };

  notifyRemoval(a);
  notifyRemoval(b);

  /* All other cases, where there may be BondStereopermutators or
   * AtomStereopermutators on a or b, should be handled correctly by
   * propagateGraphChange_.
   */

  propagateGraphChange_();
  canonicalComponentsOption_ = boost::none;
}

bool Molecule::Impl::removePermutator(const BondIndex& bond) {
  const bool removed = stereopermutators_.remove(bond);
  propagateGraphChange_();
  return removed;
}

bool Molecule::Impl::setBondType(const AtomIndex a, const AtomIndex b, const BondType bondType) {
  if(!isValidIndex_(a) || !isValidIndex_(b)) {
    throw std::out_of_range(
      "Molecule::setBondType: A supplied index is invalid!");
  }

  if(bondType == BondType::Eta) {
    throw std::logic_error("Do not manually change eta bond types, this "
                           "dynamism is handled internally");
  }

  PrivateGraph& inner = adjacencies_.inner();

  auto edgeOption = inner.edgeOption(a, b);
  if(!edgeOption) {
    addBond(a, b, bondType);
    return false;
  }

  inner.bondType(edgeOption.value()) = bondType;
  propagateGraphChange_();
  canonicalComponentsOption_ = boost::none;
  return true;
}

void Molecule::Impl::setElementType(const AtomIndex a, const Utils::ElementType elementType) {
  if(!isValidIndex_(a)) {
    throw std::out_of_range("Molecule::setElementType: This index is invalid!");
  }

  adjacencies_.inner().elementType(a) = elementType;
  propagateGraphChange_();
  canonicalComponentsOption_ = boost::none;
}

void Molecule::Impl::setShapeAtAtom(const AtomIndex a, const Shapes::Shape shape) {
  if(!isValidIndex_(a)) {
    throw std::out_of_range(
      "Molecule::setShapeAtAtom: Supplied atom index is invalid");
  }

  auto stereopermutatorOption = stereopermutators_.option(a);

  // If there is no stereopermutator at this position yet, we have to create it
  if(!stereopermutatorOption) {
    RankingInformation localRanking = rankPriority(a);

    if(localRanking.sites.size() != Shapes::size(shape)) {
      throw std::logic_error(
        "Molecule::setShapeAtAtom: The size of the supplied shape is not "
        " the same as the number of determined sites");
    }

    // Add the stereopermutator irrespective of how many assignments it has
    AtomStereopermutator newStereopermutator{
      a,
      shape,
      std::move(localRanking),
      Stereopermutators::Feasible::Functor(adjacencies_),
      AtomStereopermutator::thermalizationFunctor(adjacencies_)};

    // Default-assign stereopermutators with only one assignment
    if(newStereopermutator.numAssignments() == 1) {
      newStereopermutator.assign(0U);
    }

    stereopermutators_.add(std::move(newStereopermutator));

    propagateGraphChange_();
    canonicalComponentsOption_ = boost::none;
    return;
  }

  if(Shapes::size(stereopermutatorOption->getShape()) != Shapes::size(shape)) {
    throw std::logic_error(
      "Molecule::setShapeAtAtom: The size of the supplied shape is "
      "not the same as that of the existing stereopermutator's current shape!");
  }

  if(stereopermutatorOption->getShape() == shape) {
    // Do nothing
    return;
  }

  stereopermutatorOption->setShape(
    shape,
    Stereopermutators::Feasible::Functor(adjacencies_),
    AtomStereopermutator::thermalizationFunctor(adjacencies_));
  if(stereopermutatorOption->numAssignments() == 1) {
    stereopermutatorOption->assign(0);
  }

  // Remove any adjacent bond stereopermutators since there is no propagation
  for(BondIndex bond : adjacencies_.bonds(a)) {
    stereopermutators_.remove(bond);
  }

  propagateGraphChange_();
  canonicalComponentsOption_ = boost::none;
}

void Molecule::Impl::thermalizeStereopermutator(const AtomIndex a, const bool thermalization) {
  if(!isValidIndex_(a)) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  auto stereopermutatorOption = stereopermutators_.option(a);

  if(!stereopermutatorOption) {
    throw std::out_of_range("No stereopermutator at this index!");
  }

  stereopermutatorOption->thermalize(thermalization);

  if(thermalization) {
    /* Remove any adjacent bond stereopermutators since both constituting atom
     * stereopermutators must be assigned
     */
    for(BondIndex bond : adjacencies_.bonds(a)) {
      stereopermutators_.remove(bond);
    }
  }

  propagateGraphChange_();
  canonicalComponentsOption_ = boost::none;
}

/* Information */
boost::optional<AtomEnvironmentComponents> Molecule::Impl::canonicalComponents() const {
  return canonicalComponentsOption_;
}

std::string Molecule::Impl::str() const {
  std::stringstream info;

  if(!stereopermutators_.empty()) {
    info << "Stereopermutator information:\n";

    for(const auto& stereopermutator : stereopermutators_.atomStereopermutators()) {
      info << stereopermutator.info() << "\n";
    }

    for(const auto& stereopermutator : stereopermutators_.bondStereopermutators()) {
      info << stereopermutator.info() << "\n";
    }
  }

  return info.str();
}

boost::optional<Shapes::Shape>
Molecule::Impl::inferShape(const AtomIndex index, const RankingInformation& ranking) const {
  if(!isValidIndex_(index)) {
    throw std::out_of_range("Molecule::inferShape: Supplied index is invalid!");
  }

  if(graph().degree(index) <= 1) {
    throw std::logic_error(
      "Molecule::inferShape: No geometries exist for terminal atoms");
  }

  return ShapeInference::inferShape(adjacencies_, index, ranking);
}

std::string Molecule::Impl::dumpGraphviz() const {
  MolGraphWriter propertyWriter(&adjacencies_.inner(), &stereopermutators_);

  std::stringstream graphvizStream;

  boost::write_graphviz(
    graphvizStream, adjacencies_.inner().bgl(), propertyWriter, propertyWriter, propertyWriter);

  return graphvizStream.str();
}

const Graph& Molecule::Impl::graph() const {
  return adjacencies_;
}

std::size_t Molecule::Impl::hash() const {
  if(canonicalComponentsOption_ == boost::none) {
    throw std::logic_error("Trying to hash a non-canonical molecule.");
  }

  auto hashes = Hashes::generate(
    graph().inner(), stereopermutators_, canonicalComponentsOption_.value());

  // Convolute all of the wide hashes into a size_t hash
  static_assert(
    std::is_same<Hashes::WideHashType, boost::multiprecision::uint128_t>::value,
    "WideHash is no longer the boost multiprecision 128 uint");
  constexpr unsigned wideHashBytes = 128 / 8;
  std::vector<std::size_t> wideHashParts(wideHashBytes / sizeof(std::size_t));
  std::size_t hash = 0;
  for(const auto& wideHash : hashes) {
    boost::multiprecision::export_bits(
      wideHash, std::begin(wideHashParts), 8 * sizeof(std::size_t));
    for(const std::size_t& v : wideHashParts) {
      boost::hash_combine(hash, v);
    }
  }
  return hash;
}

const StereopermutatorList& Molecule::Impl::stereopermutators() const {
  return stereopermutators_;
}

StereopermutatorList Molecule::Impl::inferStereopermutatorsFromPositions(
  const AngstromPositions& wrapper,
  const boost::optional<std::vector<BondIndex>>& explicitBondStereopermutatorCandidatesOption,
  boost::optional<const std::unordered_set<AtomIndex>&> uninteresting,
  const boost::optional<SubstitutionsGenerator::SubstitutionMap>& substitutions) const {
  StereopermutatorList stereopermutators;

  for(const AtomIndex vertex : graph().atoms()) {
    if(uninteresting && uninteresting->count(vertex) > 0) {
      continue;
    }

    if(auto permutator = makePermutator(vertex, stereopermutators, wrapper, substitutions)) {
      stereopermutators.add(std::move(permutator.value()));
    }
  }

  const auto tryInstantiateBondStereopermutator = [&](const BondIndex& bond) {
    // Do not place stereopermutators on eta bonds.
    if(graph().bondType(bond) == BondType::Eta)
      return ;
    auto permutator = makePermutator(bond, stereopermutators, wrapper, substitutions);

    if(permutator && permutator->assigned() != boost::none) {
      stereopermutators.add(std::move(permutator.value()));
    }
  };

  if(explicitBondStereopermutatorCandidatesOption) {
    for(const BondIndex& bondIndex : *explicitBondStereopermutatorCandidatesOption) {
      tryInstantiateBondStereopermutator(bondIndex);
    }
  } else {
    // Every multiple-order bond type is a candidate
    for(const BondIndex& bondIndex : graph().bonds()) {
      if(isGraphBasedBondStereopermutatorCandidate_(graph().bondType(bondIndex))) {
        tryInstantiateBondStereopermutator(bondIndex);
      }
    }
  }

  // Look through all cycles of the molecule and try to find flat cycles
  for(const auto& cycleBonds : graph().cycles()) {

    if(cycleBonds.size() == 3) {
      continue;
    }

    const auto& cycleIndices = makeRingIndexSequence(cycleBonds);
    const double rmsPlaneDeviation
      = Cartesian::planeOfBestFitRmsd(wrapper.positions, cycleIndices);

    // Threshold for planarity
    constexpr double flatRmsPlaneDeviation = 0.07;
    if(rmsPlaneDeviation <= flatRmsPlaneDeviation) {
      for(const BondIndex& bond : cycleBonds) {
        if(!stereopermutators.option(bond)) {
          tryInstantiateBondStereopermutator(bond);
        }
      }
    }
  }
  return stereopermutators;
}

bool Molecule::Impl::canonicalCompare(
  const Molecule::Impl& other,
  const AtomEnvironmentComponents componentBitmask) const {
  /* Make sure that the components used to canonicalize each molecule instance
   * are enough to compare them canonically, too
   */
  if(
    canonicalComponentsOption_ && other.canonicalComponentsOption_
    && (canonicalComponentsOption_.value() < componentBitmask
        || other.canonicalComponentsOption_.value() < componentBitmask)) {
    throw std::logic_error("Fewer components were used in canonicalizing a "
                           "Molecule than are being compared!");
  }

  if(graph().V() != other.graph().V() || graph().E() != other.graph().E()) {
    return false;
  }

  return (
    Hashes::identityCompare(
      graph().inner(), stereopermutators(), other.graph().inner(), other.stereopermutators(), componentBitmask)
    && graph().inner().identicalGraph(other.graph().inner()));
}

boost::optional<std::vector<AtomIndex>> Molecule::Impl::modularIsomorphism(
  const Molecule::Impl& other,
  const AtomEnvironmentComponents componentBitmask) const {
  const unsigned thisNumAtoms = graph().V();

  if(thisNumAtoms != other.graph().V()) {
    return boost::none;
  }

  // Compare number of bonds too
  if(graph().E() != other.graph().E()) {
    return boost::none;
  }

  // Shortcut with same graph test if componentBitmask matches canonical
  // components
  if(canonicalComponents() == componentBitmask && other.canonicalComponents() == componentBitmask) {
    if(canonicalCompare(other, componentBitmask)) {
      return Temple::iota<AtomIndex>(thisNumAtoms);
    }

    return boost::none;
  }

  /* boost isomorphism will allocate a vector of size maxHash, this is dangerous
   * as the maximum hash can be immense, another post-processing step is needed
   * for the calculated hashes to decrease the memory space requirements
   *
   * This maps the hashes to an incremented number:
   */
  std::vector<Hashes::HashType> thisHashes;
  std::vector<Hashes::HashType> otherHashes;
  Hashes::HashType maxHash = 0;

  std::tie(thisHashes, otherHashes, maxHash) = Hashes::narrow(
    Hashes::generate(graph().inner(), stereopermutators(), componentBitmask),
    Hashes::generate(other.graph().inner(), other.stereopermutators(), componentBitmask));

  // Where the corresponding index from the other graph is stored
  std::vector<AtomIndex> indexMap(thisNumAtoms);

  const auto& thisBGL = adjacencies_.inner().bgl();
  const auto& otherBGL = other.adjacencies_.inner().bgl();

  const bool isomorphic = boost::isomorphism(
    thisBGL,
    otherBGL,
    boost::make_safe_iterator_property_map(
      indexMap.begin(), thisNumAtoms, boost::get(boost::vertex_index, thisBGL)),
    Hashes::LookupFunctor(thisHashes),
    Hashes::LookupFunctor(otherHashes),
    maxHash,
    boost::get(boost::vertex_index, thisBGL),
    boost::get(boost::vertex_index, otherBGL));

  if(isomorphic) {
    return indexMap;
  }

  return boost::none;
}

RankingInformation Molecule::Impl::rankPriority(
  const AtomIndex a,
  const std::vector<AtomIndex>& excludeAdjacent,
  const boost::optional<AngstromPositions>& positionsOption) const {
  if(!isValidIndex_(a)) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  RankingInformation rankingResult;

  // Expects that bond types are set properly, complains otherwise
  rankingResult.sites
    = GraphAlgorithms::sites(adjacencies_.inner(), a, excludeAdjacent);

  std::string molGraphviz;
#ifndef NDEBUG
  molGraphviz = dumpGraphviz();
#endif

  // Rank the substituents
  auto expandedTree = RankingTree(
    graph(),
    stereopermutators(),
    molGraphviz,
    a,
    excludeAdjacent,
    RankingTree::ExpansionOption::OnlyRequiredBranches,
    positionsOption);

  rankingResult.substituentRanking = expandedTree.getRanked();

  // Combine site information and substituent ranking into a site ranking
  rankingResult.siteRanking = RankingInformation::rankSites(
    rankingResult.sites, rankingResult.substituentRanking);

  // Find links between sites
  rankingResult.links = GraphAlgorithms::siteLinks(
    adjacencies_.inner(), a, rankingResult.sites, excludeAdjacent);

  return rankingResult;
}

bool Molecule::Impl::operator==(const Impl& other) const {
  if(canonicalComponentsOption_ == AtomEnvironmentComponents::All && other.canonicalComponentsOption_ == AtomEnvironmentComponents::All) {
    return canonicalCompare(other, AtomEnvironmentComponents::All);
  }

  // Better with newer boost: .has_value()
  return static_cast<bool>(modularIsomorphism(other, AtomEnvironmentComponents::All));
}

bool Molecule::Impl::operator!=(const Impl& other) const {
  return !(*this == other);
}

} // namespace Molassembler
} // namespace Scine
