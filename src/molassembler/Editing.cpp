/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Editing.h"

#include "molassembler/Graph/Bridge.h"
#include "molassembler/Graph/PrivateGraph.h"
#include "molassembler/Molecule/MoleculeImpl.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/BondStereopermutator.h"
#include "molassembler/Stereopermutators/AbstractPermutations.h"
#include "molassembler/Stereopermutators/FeasiblePermutations.h"

#include "Utils/Geometry/ElementInfo.h"

#include "molassembler/Temple/Functional.h"

namespace Scine {
namespace molassembler {
namespace {

/**
 * @brief Copy a subset of vertices and edges into another graph
 *
 * @param source The source graph from which to copy
 * @param target The target graph into which to copy
 * @param copyVertices List of vertices to copy. If empty, copies all vertices.
 *
 * @return An unordered map of source atom indices to target atom indices
 */
std::unordered_map<AtomIndex, AtomIndex> transferGraph(
  const PrivateGraph& source,
  PrivateGraph& target,
  const std::vector<AtomIndex>& copyVertices
) {
  /* Copy over element types, creating the new vertices in target, collecting
   * the new vertex indices in the target
   */
  std::unordered_map<AtomIndex, AtomIndex> copyVertexTargetIndices;
  if(copyVertices.empty()) {
    for(const AtomIndex& vertex : source.vertices()) {
      copyVertexTargetIndices.insert({
        vertex,
        target.addVertex(
          source.elementType(vertex)
        )
      });
    }
  } else {
    for(const AtomIndex& vertexIndex : copyVertices) {
      copyVertexTargetIndices.insert({
        vertexIndex,
        target.addVertex(
          source.elementType(vertexIndex)
        )
      });
    }
  }

  for(const PrivateGraph::Edge& e : source.edges()) {
    auto findSourceIter = copyVertexTargetIndices.find(
      source.source(e)
    );

    auto findTargetIter = copyVertexTargetIndices.find(
      source.target(e)
    );

    /* If both source and target indices are keys in copyVertexTargetIndices,
     * then we copy the edge
     */
    if(
      findSourceIter != std::end(copyVertexTargetIndices)
      && findTargetIter != std::end(copyVertexTargetIndices)
    ) {
      target.addEdge(
        findSourceIter->second,
        findTargetIter->second,
        source.bondType(e)
      );
    }
  }

  return copyVertexTargetIndices;
}

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

} // namespace

std::pair<Molecule, Molecule> Editing::cleave(const Molecule& a, const BondIndex bridge) {
  if(a.graph().canRemove(bridge)) {
    throw std::logic_error(
      "The supplied bond can be removed without disconnecting the graph. It is "
      "not a suitable edge to cleave a molecule in two."
    );
  }

  const AtomIndex N = a.graph().N();

  // Discover which vertices belong to which component after cleaving
  auto sides = a.graph().inner().splitAlongBridge(
    toInner(bridge, a.graph().inner())
  );

  // Ensure both sides will fulfill preconditions to be a molecule (at least one atom)
  assert(!sides.first.empty() && !sides.second.empty());

  // Construct separate OuterGraphs for each component of the disconnected graph
  std::pair<PrivateGraph, PrivateGraph> graphs;
  auto vertexMappings = std::make_pair(
    transferGraph(
      a.graph().inner(),
      graphs.first,
      sides.first
    ),
    transferGraph(
      a.graph().inner(),
      graphs.second,
      sides.second
    )
  );

  /* Copy stereopermutators, adapting internal state according to index
   * mappings to either component
   */
  std::pair<StereopermutatorList, StereopermutatorList> stereopermutatorLists;
  transferStereopermutators(
    a.stereopermutators(),
    stereopermutatorLists.first,
    vertexMappings.first,
    N
  );
  transferStereopermutators(
    a.stereopermutators(),
    stereopermutatorLists.second,
    vertexMappings.second,
    N
  );

  // Make molecules out of the components
  auto molecules = std::make_pair<Molecule, Molecule>(
    Molecule(
      Graph(std::move(graphs.first)),
      std::move(stereopermutatorLists.first)
    ),
    Molecule(
      Graph(std::move(graphs.second)),
      std::move(stereopermutatorLists.second)
    )
  );

  /* Follow a procedure similar to Molecule::Impl::removeBond to get valid
   * state after removal:
   */
  auto fixEdgeAndPropagate = [](
    Molecule& molecule,
    const AtomIndex notifyIndex
  ) {
    if(auto stereopermutatorOption = molecule.pImpl_->stereopermutators_.option(notifyIndex)) {
      auto localRanking = molecule.rankPriority(notifyIndex);

      // In case the central atom becomes terminal, just drop the stereopermutator
      if(localRanking.sites.size() <= 1) {
        molecule.pImpl_->stereopermutators_.remove(notifyIndex);
        return;
      }

      boost::optional<shapes::Shape> shapeOption;
      if(Options::shapeTransition == ShapeTransition::PrioritizeInferenceFromGraph) {
        shapeOption = molecule.inferShape(notifyIndex, localRanking);
      }

      // Notify the stereopermutator to remove the placeholder
      stereopermutatorOption->propagate(
        molecule.pImpl_->adjacencies_,
        std::move(localRanking),
        shapeOption
      );

      // Default-assign if possible
      if(
        stereopermutatorOption->assigned() == boost::none
        && stereopermutatorOption->numStereopermutations() == 1
        && stereopermutatorOption->numAssignments() == 1
      ) {
        stereopermutatorOption->assign(0);
      }

      // TODO notify or just remove any bond stereopermutators on notifyIndex
    }

    // Rerank everywhere and update all stereopermutators
    molecule.pImpl_->propagateGraphChange_();
  };

  fixEdgeAndPropagate(
    molecules.first,
    vertexMappings.first.at(bridge.first)
  );

  fixEdgeAndPropagate(
    molecules.second,
    vertexMappings.second.at(bridge.second)
  );

  return molecules;
}

Molecule Editing::insert(
  Molecule log,
  const Molecule& wedge,
  const BondIndex logBond,
  const AtomIndex firstWedgeAtom,
  const AtomIndex secondWedgeAtom
) {
  const AtomIndex logN = log.graph().N();

  /* - Batch insert wedge appropriately
   * - Disconnect log at logBond and connect wedge atoms
   * - Copy in stereopermutators from wedge with an index mapping adjusting
   *   internal state
   * - Call notifyGraphChange and return log
   */
  PrivateGraph& logInner = log.pImpl_->adjacencies_.inner();

  // Copy all vertices from wedge into log
  auto vertexMapping = transferGraph(
    wedge.graph().inner(),
    logInner,
    {}
  );

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
  StereopermutatorList& logStereopermutators = log.pImpl_->stereopermutators_;

  transferStereopermutators(
    wedge.stereopermutators(),
    logStereopermutators,
    vertexMapping,
    wedge.graph().N()
  );

  /* Two things remain to be done for each new bond:
   * - Log's stereopermutator ranking needs to be updated with the new atom
   *   index of its substituent
   * - Wedge's stereopermutator on the new index for the wedge index needs to
   *   be notified that it has a new substituent
   */
  auto iota = temple::iota<AtomIndex>(logN);
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

      boost::optional<shapes::Shape> shapeOption;
      if(Options::shapeTransition == ShapeTransition::PrioritizeInferenceFromGraph) {
        shapeOption = log.inferShape(newWedgeIndex, localRanking);
      }

      permutatorOption->propagate(
        log.pImpl_->adjacencies_,
        std::move(localRanking),
        shapeOption
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
  log.pImpl_->propagateGraphChange_();
  return log;
}

Molecule Editing::superpose(
  Molecule top,
  const Molecule& bottom,
  const AtomIndex topAtom,
  const AtomIndex bottomAtom
) {
  // Copy in all vertices except bottomAtom from bottom into top
  std::vector<AtomIndex> bottomCopyAtoms(bottom.graph().N() - 1);
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

  PrivateGraph& topInner = top.pImpl_->adjacencies_.inner();

  auto vertexMapping = transferGraph(
    bottom.graph().inner(),
    topInner,
    bottomCopyAtoms
  );

  /* Now for some trickery:
   * We can use transferStereopermutators to adjust the internal state of
   * permutators on and adjacent to bottomAtom during the copy if we change
   * the vertex mapping to map bottomAtom to topAtom. Then we don't have to
   * manually adjust them later.
   */

  // Copy in stereopermutators from bottom, including ones placed on bottomAtom
  StereopermutatorList& topStereopermutators = top.pImpl_->stereopermutators_;
  vertexMapping[bottomAtom] = topAtom;
  transferStereopermutators(
    bottom.stereopermutators(),
    topStereopermutators,
    vertexMapping,
    bottom.graph().N(),
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

      boost::optional<shapes::Shape> shapeOption;
      if(Options::shapeTransition == ShapeTransition::PrioritizeInferenceFromGraph) {
        shapeOption = top.inferShape(topAtom, localRanking);
      }

      topPermutatorOption->propagate(
        top.pImpl_->adjacencies_,
        std::move(localRanking),
        shapeOption
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
  top.pImpl_->propagateGraphChange_();
  return top;
}

Molecule Editing::substitute(
  const Molecule& left,
  const Molecule& right,
  const BondIndex leftBond,
  const BondIndex rightBond
) {
  PrivateGraph innerGraph;
  StereopermutatorList stereopermutators;

  // Identify sides of each bond
  auto leftSides = left.graph().splitAlongBridge(leftBond);
  auto rightSides = right.graph().splitAlongBridge(rightBond);

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

    double aWeight = temple::accumulate(
      sideA,
      0.0,
      [&molecule](double carry, const AtomIndex i) -> double {
        return carry + Utils::ElementInfo::mass(molecule.graph().elementType(i));
      }
    );

    double bWeight = temple::accumulate(
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
  AtomIndex leftHeavierBondSide, leftLighterBondSide, rightHeavierBondSide, rightLighterBondSide;

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

  assert(
    temple::makeContainsPredicate(leftHeavierSide)(leftHeavierBondSide)
  );
  assert(
    temple::makeContainsPredicate(rightHeavierSide)(rightHeavierBondSide)
  );

  // Copy over graphs
  auto leftVertexMapping = transferGraph(
    left.graph().inner(),
    innerGraph,
    leftHeavierSide
  );

  auto rightVertexMapping = transferGraph(
    right.graph().inner(),
    innerGraph,
    rightHeavierSide
  );

  // Copy over left stereopermutators
  leftVertexMapping[leftLighterBondSide] = rightVertexMapping.at(rightHeavierBondSide);
  transferStereopermutators(
    left.stereopermutators(),
    stereopermutators,
    leftVertexMapping,
    left.graph().N(),
    {leftLighterBondSide}
  );

  // Copy over right stereopermutators
  rightVertexMapping[rightLighterBondSide] = leftVertexMapping.at(leftHeavierBondSide);
  transferStereopermutators(
    right.stereopermutators(),
    stereopermutators,
    rightVertexMapping,
    right.graph().N(),
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
  compound.pImpl_->propagateGraphChange_();
  return compound;
}

Molecule Editing::connect(
  Molecule a,
  const Molecule& b,
  const AtomIndex aConnectAtom,
  const AtomIndex bConnectAtom,
  const BondType bondType
) {
  PrivateGraph& aInnerGraph = a.pImpl_->adjacencies_.inner();
  StereopermutatorList& aStereopermutators = a.pImpl_->stereopermutators_;

  // Copy b's graph into a
  auto vertexMapping = transferGraph(
    b.graph().inner(),
    aInnerGraph,
    {}
  );

  // Copy b's stereopermutators into a
  transferStereopermutators(
    b.stereopermutators(),
    aStereopermutators,
    vertexMapping,
    b.graph().N()
  );

  // Add the bond (propagating stereopermutator state and reranking everywhere)
  a.addBond(aConnectAtom, vertexMapping.at(bConnectAtom), bondType);

  return a;
}

Molecule Editing::addLigand(
  Molecule a,
  const Molecule& ligand,
  AtomIndex complexatingAtom,
  const std::vector<AtomIndex>& ligandBindingAtoms
) {
  PrivateGraph& aInnerGraph = a.pImpl_->adjacencies_.inner();
  StereopermutatorList& aStereopermutators = a.pImpl_->stereopermutators_;

  auto vertexMapping = transferGraph(
    ligand.graph().inner(),
    aInnerGraph,
    {}
  );

  // Copy b's stereopermutators into a
  transferStereopermutators(
    ligand.stereopermutators(),
    aStereopermutators,
    vertexMapping,
    ligand.graph().N()
  );

  for(const AtomIndex bindingAtom : ligandBindingAtoms) {
    a.addBond(complexatingAtom, vertexMapping.at(bindingAtom), BondType::Single);
  }

  return a;
}

} // namespace molassembler
} // namespace Scine
