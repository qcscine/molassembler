/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Editing.h"

#include "Molassembler/Graph/Bridge.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/Molecule/MoleculeImpl.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/Stereopermutators/AbstractPermutations.h"
#include "Molassembler/Stereopermutators/FeasiblePermutations.h"

#include "Utils/Geometry/ElementInfo.h"

#include "Molassembler/Temple/Functional.h"

namespace Scine {
namespace Molassembler {
namespace {

/**
 * @brief Transfer stereopermutators between stereopermutator lists
 *
 * @param source Source stereopermutator list to copy out of
 * @param target Target stereopermutator list to copy into
 * @param vertexMapping Mapping of vertex indices from the source to the target
 *   index.
 * @param sourceN Number of vertices in the molecule of which @p source is part:
 *   This is necessary to construct the index permutation that is applied to
 *   each stereopermutator. All vertices that are not contained in vertexMapping
 *   are mapped to PrivateGraph::removalPlaceholder.
 * @param doNotCopyVertices Vertex indices that are part of vertexMapping but
 *   of which stereopermutators should not be copied.
 *
 * @note If there are stereopermutators on an edge between existing
 * stereopermutators and new ones, you can propagate the internal state of the
 * new stereopermutators by populating vertexMapping with mappings from
 * non-copied atoms to atoms from existing vertices in the molecule modeled by
 * @p target and specifying for which of these no stereopermutators should be
 * copied in @p doNotCopyVertices.
 */
void transferStereopermutators(
  const StereopermutatorList& source,
  StereopermutatorList& target,
  const std::unordered_map<AtomIndex, AtomIndex>& vertexMapping,
  const AtomIndex sourceN,
  const std::unordered_set<AtomIndex>& doNotCopyVertices = {}
) {
  /* Copy over stereopermutators */
  /* How to deal with source stereopermutators on the edge? They have atom
   * indices in their state that do not exist in the target.
   *
   * Can probably reuse their applyPermutation methods if I can figure out the
   * right permutation.
   */

  /* Construct a permutation that resolves atom indices from sources using the
   * vertex mapping
   */
  std::vector<AtomIndex> permutation (sourceN, PrivateGraph::removalPlaceholder);
  for(const auto& iterPair : vertexMapping) {
    permutation.at(iterPair.first) = iterPair.second;
  }

  // Copy over permutators, adjusting their internal state
  for(const auto& permutator : source.atomStereopermutators()) {
    if(
      vertexMapping.count(permutator.placement()) > 0
      && doNotCopyVertices.count(permutator.placement()) == 0
    ) {
      AtomStereopermutator copy = permutator;
      copy.applyPermutation(permutation);
      target.add(std::move(copy));
    }
  }

  for(const auto& permutator : source.bondStereopermutators()) {
    BondIndex edge = permutator.placement();
    if(
      vertexMapping.count(edge.first) > 0
      && doNotCopyVertices.count(edge.first) == 0
      && vertexMapping.count(edge.second) > 0
      && doNotCopyVertices.count(edge.second) == 0
    ) {
      BondStereopermutator copy = permutator;
      copy.applyPermutation(permutation);
      target.add(std::move(copy));
    }
  }
}

Editing::Cleaved cleaveImpl(
  const Molecule& a,
  const AtomIndex left,
  const std::vector<AtomIndex>& right,
  const std::vector<AtomIndex>& leftSide,
  const std::vector<AtomIndex>& rightSide
) {
  const AtomIndex V = a.graph().V();

  // Construct separate OuterGraphs for each component of the disconnected graph
  std::pair<PrivateGraph, PrivateGraph> graphs;
  const auto vertexMappings = std::make_pair(
    graphs.first.merge(a.graph().inner(), leftSide),
    graphs.second.merge(a.graph().inner(), rightSide)
  );

  /* Copy stereopermutators, adapting internal state according to index
   * mappings to either component
   */
  std::pair<StereopermutatorList, StereopermutatorList> stereopermutatorLists;
  transferStereopermutators(
    a.stereopermutators(),
    stereopermutatorLists.first,
    vertexMappings.first,
    V
  );
  transferStereopermutators(
    a.stereopermutators(),
    stereopermutatorLists.second,
    vertexMappings.second,
    V
  );

  // Make molecules out of the components
  Editing::Cleaved cleaved {
    Molecule(
      Graph(std::move(graphs.first)),
      std::move(stereopermutatorLists.first)
    ),
    Molecule(
      Graph(std::move(graphs.second)),
      std::move(stereopermutatorLists.second)
    ),
    {}
  };

  // Write the component map
  cleaved.componentMap.resize(V);
  for(const auto& firstComponentPair : vertexMappings.first) {
    cleaved.componentMap.at(firstComponentPair.first) = std::make_pair(
      0,
      firstComponentPair.second
    );
  }
  for(const auto& secondComponentPair : vertexMappings.second) {
    cleaved.componentMap.at(secondComponentPair.first) = std::make_pair(
      1,
      secondComponentPair.second
    );
  }

  /* Follow a procedure similar to Molecule::Impl::removeBond to get valid
   * state after removal:
   */
  auto fixEdgeAndPropagate = [](
    Molecule& molecule,
    const AtomIndex notifyIndex
  ) {
    if(auto stereopermutatorOption = molecule.stereopermutators(Molecule::unsafe_tag).option(notifyIndex)) {
      auto localRanking = molecule.rankPriority(notifyIndex);

      // In case the central atom becomes terminal, just drop the stereopermutator
      if(localRanking.sites.size() <= 1) {
        molecule.stereopermutators(Molecule::unsafe_tag).remove(notifyIndex);
        return;
      }

      boost::optional<Shapes::Shape> shapeOption;
      if(Options::shapeTransition == ShapeTransition::PrioritizeInferenceFromGraph) {
        shapeOption = molecule.inferShape(notifyIndex, localRanking);
      }

      // Notify the stereopermutator to remove the placeholder
      stereopermutatorOption->propagate(
        std::move(localRanking),
        shapeOption,
        Stereopermutators::Feasible::Functor(molecule.graph()),
        AtomStereopermutator::thermalizationFunctor(molecule.graph())
      );

      // Default-assign if possible
      if(
        stereopermutatorOption->assigned() == boost::none
        && stereopermutatorOption->numStereopermutations() == 1
        && stereopermutatorOption->numAssignments() == 1
      ) {
        stereopermutatorOption->assign(0);
      }

      /* TODO propagate bond stereopermutators on the adjacent vertex instead
       * of dropping them
       */
      for(const BondIndex bond : molecule.graph().bonds(notifyIndex)) {
        if(molecule.stereopermutators().option(bond)) {
          molecule.stereopermutators(Molecule::unsafe_tag).remove(bond);
        }
      }
    }

    // Rerank everywhere and update all stereopermutators
    molecule.propagate(Molecule::unsafe_tag);
  };

  for(const AtomIndex i : right) {
    if(a.graph().bondType(BondIndex {left, i}) != BondType::Eta) {
      fixEdgeAndPropagate(
        cleaved.second,
        vertexMappings.second.at(i)
      );
    }
  }

  fixEdgeAndPropagate(
    cleaved.first,
    vertexMappings.first.at(left)
  );

  return cleaved;
}

} // namespace

namespace Editing {

Cleaved cleave(
  const Molecule& a,
  const AtomSitePair site
) {
  auto maybePermutator = a.stereopermutators().option(site.first);
  if(!maybePermutator) {
    throw std::invalid_argument("No stereopermutator on specified atom");
  }
  const auto& siteAtoms = maybePermutator->getRanking().sites.at(site.second);
  const auto sides = a.graph().inner().splitAlongBridge(
    site.first,
    siteAtoms
  );
  assert(!sides.first.empty() && !sides.second.empty());
  return cleaveImpl(
    a,
    site.first,
    siteAtoms,
    sides.first,
    sides.second
  );
}

Cleaved cleave(const Molecule& a, const BondIndex bridge) {
  if(a.graph().canRemove(bridge)) {
    throw std::logic_error(
      "The supplied bond can be removed without disconnecting the graph. It is "
      "not a suitable edge to cleave a molecule in two."
    );
  }

  // Discover which vertices belong to which component after cleaving
  const auto sides = a.graph().inner().splitAlongBridge(
    toInner(bridge, a.graph().inner())
  );

  assert(!sides.first.empty() && !sides.second.empty());
  Cleaved cleaved = cleaveImpl(
    a,
    bridge.first,
    std::vector<AtomIndex>(1, bridge.second),
    sides.first,
    sides.second
  );
  return cleaved;
}

Molecule insert(
  Molecule log,
  const Molecule& wedge,
  const BondIndex logBond,
  const AtomIndex firstWedgeAtom,
  const AtomIndex secondWedgeAtom
) {
  const AtomIndex logN = log.graph().V();

  /* - Batch insert wedge appropriately
   * - Disconnect log at logBond and connect wedge atoms
   * - Copy in stereopermutators from wedge with an index mapping adjusting
   *   internal state
   * - Call notifyGraphChange and return log
   */
  PrivateGraph& logInner = log.graph(Molecule::unsafe_tag).inner();

  // Copy all vertices from wedge into log
  auto vertexMapping = logInner.merge(wedge.graph().inner());

  // Disconnect log at logBond
  PrivateGraph::Edge logEdge = toInner(logBond, logInner);
  BondType bondType = logInner.bondType(logEdge);
  logInner.removeEdge(logEdge);

  // Connect the wedge atoms
  logInner.addEdge(
    logInner.source(logEdge),
    vertexMapping.at(firstWedgeAtom),
    bondType
  );

  logInner.addEdge(
    logInner.target(logEdge),
    vertexMapping.at(secondWedgeAtom),
    bondType
  );

  // Copy in stereopermutators from wedge
  StereopermutatorList& logStereopermutators = log.stereopermutators(Molecule::unsafe_tag);

  transferStereopermutators(
    wedge.stereopermutators(),
    logStereopermutators,
    vertexMapping,
    wedge.graph().V()
  );

  /* Two things remain to be done for each new bond:
   * - Log's stereopermutator ranking needs to be updated with the new atom
   *   index of its substituent
   * - Wedge's stereopermutator on the new index for the wedge index needs to
   *   be notified that it has a new substituent
   */
  auto iota = Temple::iota<AtomIndex>(logN);
  auto fixNewBond = [&](
    AtomIndex logSide,
    AtomIndex otherLogSide,
    AtomIndex oldWedgeIndex
  ) {
    const AtomIndex newWedgeIndex = vertexMapping.at(oldWedgeIndex);
    /* Modify iota so that it becomes a permutation that merely changes the
     * old bonded atom index into the new one
     */
    auto& permutation = iota;
    permutation.at(otherLogSide) = newWedgeIndex;

    auto stereopermutatorOption = logStereopermutators.option(logSide);
    assert(stereopermutatorOption);

    stereopermutatorOption->applyPermutation(permutation);

    // Re-establish permutation as iota
    permutation.at(otherLogSide) = otherLogSide;

    // Notify the new wedge stereopermutator on the other end of the bond to logSide
    if(auto permutatorOption = logStereopermutators.option(newWedgeIndex)) {
      auto localRanking = log.rankPriority(newWedgeIndex);

      boost::optional<Shapes::Shape> shapeOption;
      if(Options::shapeTransition == ShapeTransition::PrioritizeInferenceFromGraph) {
        shapeOption = log.inferShape(newWedgeIndex, localRanking);
      }

      permutatorOption->propagate(
        std::move(localRanking),
        shapeOption,
        Stereopermutators::Feasible::Functor(log.graph()),
        AtomStereopermutator::thermalizationFunctor(log.graph())
      );

      // Default assign if possible
      if(
        permutatorOption->assigned() == boost::none
        && permutatorOption->numStereopermutations() == 1
        && permutatorOption->numAssignments() == 1
      ) {
        permutatorOption->assign(0);
      }

      // TODO notify or just remove bond stereopermutators on the same index
    }
  };

  fixNewBond(
    logBond.first,
    logBond.second,
    firstWedgeAtom
  );

  fixNewBond(
    logBond.second,
    logBond.first,
    secondWedgeAtom
  );

  // Rerank everywhere and return
  log.propagate(Molecule::unsafe_tag);
  return log;
}

Molecule superpose(
  Molecule top,
  const Molecule& bottom,
  const AtomIndex topAtom,
  const AtomIndex bottomAtom
) {
  // Copy in all vertices except bottomAtom from bottom into top
  std::vector<AtomIndex> bottomCopyAtoms(bottom.graph().V() - 1);
  std::iota(
    std::begin(bottomCopyAtoms),
    std::begin(bottomCopyAtoms) + bottomAtom,
    0
  );
  std::iota(
    std::begin(bottomCopyAtoms) + bottomAtom,
    std::end(bottomCopyAtoms),
    bottomAtom + 1
  );

  PrivateGraph& topInner = top.graph(Molecule::unsafe_tag).inner();

  auto vertexMapping = topInner.merge(
    bottom.graph().inner(),
    bottomCopyAtoms
  );

  /* Now for some trickery:
   * We can use transferStereopermutators to adjust the internal state of
   * permutators on and adjacent to bottomAtom during the copy if we change
   * the vertex mapping to map bottomAtom to topAtom. Then we don't have to
   * manually adjust them later.
   */

  // Copy in stereopermutators from bottom, including ones placed on bottomAtom
  StereopermutatorList& topStereopermutators = top.stereopermutators(Molecule::unsafe_tag);
  vertexMapping[bottomAtom] = topAtom;
  transferStereopermutators(
    bottom.stereopermutators(),
    topStereopermutators,
    vertexMapping,
    bottom.graph().V(),
    {bottomAtom}
  );

  /* Copy in adjacencies of bottomAtom manually and notify the stereopermutator
   * on top each time
   */
  const PrivateGraph& bottomInner = bottom.graph().inner();

  auto topPermutatorOption = topStereopermutators.option(topAtom);

  for(const AtomIndex bottomAtomAdjacent : bottomInner.adjacents(bottomAtom)) {
    const AtomIndex newSubstituent = vertexMapping.at(bottomAtomAdjacent);
    const BondType bondType = bottomInner.bondType(
      bottomInner.edge(bottomAtom, bottomAtomAdjacent)
    );
    topInner.addEdge(
      topAtom,
      newSubstituent,
      bondType
    );

    if(topPermutatorOption) {
      auto localRanking = top.rankPriority(topAtom);

      boost::optional<Shapes::Shape> shapeOption;
      if(Options::shapeTransition == ShapeTransition::PrioritizeInferenceFromGraph) {
        shapeOption = top.inferShape(topAtom, localRanking);
      }

      topPermutatorOption->propagate(
        std::move(localRanking),
        shapeOption,
        Stereopermutators::Feasible::Functor(top.graph()),
        AtomStereopermutator::thermalizationFunctor(top.graph())
      );

      // Default assign if possible
      if(
        topPermutatorOption->assigned() == boost::none
        && topPermutatorOption->numStereopermutations() == 1
        && topPermutatorOption->numAssignments() == 1
      ) {
        topPermutatorOption->assign(0);
      }
    }
  }

  // Rerank everywhere and return
  top.propagate(Molecule::unsafe_tag);
  return top;
}

Molecule substitute(
  const Molecule& left,
  const Molecule& right,
  const BondIndex leftBond,
  const BondIndex rightBond
) {
  PrivateGraph innerGraph;
  StereopermutatorList stereopermutators;

  // Identify sides of each bond
  const auto leftSides = left.graph().splitAlongBridge(leftBond);
  const auto rightSides = right.graph().splitAlongBridge(rightBond);

  // Figure out which side of each bond is 'heavier'
  auto sideCompare = [](
    const Molecule& molecule,
    const std::vector<AtomIndex>& sideA,
    const std::vector<AtomIndex>& sideB
  ) -> bool {
    if(sideA.size() < sideB.size()) {
      return true;
    }

    if(sideA.size() > sideB.size()) {
      return false;
    }

    const double aWeight = Temple::accumulate(
      sideA,
      0.0,
      [&molecule](double carry, const AtomIndex i) -> double {
        return carry + Utils::ElementInfo::mass(molecule.graph().elementType(i));
      }
    );

    const double bWeight = Temple::accumulate(
      sideB,
      0.0,
      [&molecule](double carry, const AtomIndex i) -> double {
        return carry + Utils::ElementInfo::mass(molecule.graph().elementType(i));
      }
    );

    return aWeight < bWeight;
  };

  const auto& leftHeavierSide = sideCompare(left, leftSides.first, leftSides.second) ? leftSides.second : leftSides.first;
  const auto& rightHeavierSide = sideCompare(right, rightSides.first, rightSides.second) ? rightSides.second : rightSides.first;

  // Figure out which bond index of the bond is the one that belongs to the heavier side
  AtomIndex leftHeavierBondSide {};
  AtomIndex leftLighterBondSide {};
  AtomIndex rightHeavierBondSide {};
  AtomIndex rightLighterBondSide {};

  if(std::addressof(leftHeavierSide) == std::addressof(leftSides.first)) {
    leftHeavierBondSide = leftBond.first;
    leftLighterBondSide = leftBond.second;
  } else {
    leftHeavierBondSide = leftBond.second;
    leftLighterBondSide = leftBond.first;
  }

  if(std::addressof(rightHeavierSide) == std::addressof(rightSides.first)) {
    rightHeavierBondSide = rightBond.first;
    rightLighterBondSide = rightBond.second;
  } else {
    rightHeavierBondSide = rightBond.second;
    rightLighterBondSide = rightBond.first;
  }

  assert(Temple::makeContainsPredicate(leftHeavierSide)(leftHeavierBondSide));
  assert(Temple::makeContainsPredicate(rightHeavierSide)(rightHeavierBondSide));

  // Copy over graphs
  auto leftVertexMapping = innerGraph.merge(
    left.graph().inner(),
    leftHeavierSide
  );

  auto rightVertexMapping = innerGraph.merge(
    right.graph().inner(),
    rightHeavierSide
  );

  // Copy over left stereopermutators
  leftVertexMapping[leftLighterBondSide] = rightVertexMapping.at(rightHeavierBondSide);
  transferStereopermutators(
    left.stereopermutators(),
    stereopermutators,
    leftVertexMapping,
    left.graph().V(),
    {leftLighterBondSide}
  );

  // Copy over right stereopermutators
  rightVertexMapping[rightLighterBondSide] = leftVertexMapping.at(leftHeavierBondSide);
  transferStereopermutators(
    right.stereopermutators(),
    stereopermutators,
    rightVertexMapping,
    right.graph().V(),
    {rightLighterBondSide}
  );

  // Add the missing bond
  innerGraph.addEdge(
    leftVertexMapping.at(leftHeavierBondSide),
    rightVertexMapping.at(rightHeavierBondSide),
    BondType::Single
  );

  // Make a molecule out of the components
  Molecule compound {
    Graph(std::move(innerGraph)),
    std::move(stereopermutators)
  };

  // Rerank everywhere and return
  compound.propagate(Molecule::unsafe_tag);
  return compound;
}

Molecule connect(
  Molecule a,
  const Molecule& b,
  const AtomIndex aConnectAtom,
  const AtomIndex bConnectAtom,
  const BondType bondType
) {
  PrivateGraph& aInnerGraph = a.graph(Molecule::unsafe_tag).inner();
  StereopermutatorList& aStereopermutators = a.stereopermutators(Molecule::unsafe_tag);

  // Copy b's graph into a
  auto vertexMapping = aInnerGraph.merge(b.graph().inner());

  // Copy b's stereopermutators into a
  transferStereopermutators(
    b.stereopermutators(),
    aStereopermutators,
    vertexMapping,
    b.graph().V()
  );

  // Add the bond (propagating stereopermutator state and reranking everywhere)
  a.addBond(aConnectAtom, vertexMapping.at(bConnectAtom), bondType);

  return a;
}

Molecule addLigand(
  Molecule a,
  const Molecule& ligand,
  AtomIndex complexatingAtom,
  const std::vector<AtomIndex>& ligandBindingAtoms
) {
  PrivateGraph& aInnerGraph = a.graph(Molecule::unsafe_tag).inner();
  StereopermutatorList& aStereopermutators = a.stereopermutators(Molecule::unsafe_tag);

  auto vertexMapping = aInnerGraph.merge(ligand.graph().inner());

  // Copy b's stereopermutators into a
  transferStereopermutators(
    ligand.stereopermutators(),
    aStereopermutators,
    vertexMapping,
    ligand.graph().V()
  );

  for(const AtomIndex bindingAtom : ligandBindingAtoms) {
    a.addBond(complexatingAtom, vertexMapping.at(bindingAtom), BondType::Single);
  }

  return a;
}

} // namespace Editing
} // namespace Molassembler
} // namespace Scine
