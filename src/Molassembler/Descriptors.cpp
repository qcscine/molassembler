/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Descriptors.h"

#include "Molassembler/Graph.h"
#include "Molassembler/Cycles.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Temple/Adaptors/SequentialPairs.h"
#include "Molassembler/Temple/Functional.h"

#include "boost/bimap.hpp"

namespace Scine {
namespace Molassembler {
namespace {

struct DescendantSubgraph {
  using Bimap = boost::bimap<AtomIndex, AtomIndex>;

  DescendantSubgraph(
    const AtomIndex source,
    const std::vector<AtomIndex>& descendants,
    const PrivateGraph& molGraph
  ) {
    // Add the source as anchor
    const AtomIndex N = descendants.size();
    for(AtomIndex i = 0; i < N; ++i) {
      if(descendants.at(i) == source) {
        const AtomIndex newVertex = subgraph.addVertex(molGraph.elementType(i));
        indexMap.insert(Bimap::value_type(i, newVertex));
      }
    }

    // Copy edges
    for(const auto& vertexPair : indexMap.left) {
      const AtomIndex molGraphI = vertexPair.first;
      for(const AtomIndex molGraphJ : molGraph.adjacents(molGraphI)) {
        const auto findIter = indexMap.left.find(molGraphJ);
        if(findIter != indexMap.left.end()) {
          if(!subgraph.edgeOption(vertexPair.second, findIter->second)) {
            subgraph.addEdge(
              vertexPair.second,
              findIter->second,
              molGraph.bondType(molGraph.edge(molGraphI, molGraphJ))
            );
          }
        }
      }
    }

    if(subgraph.connectedComponents() != 1) {
      throw std::runtime_error("Subgraph construction of descendant did not yield a single connected component!");
    }
  }

  PrivateGraph subgraph;
  Bimap indexMap;
};

} // namespace

unsigned numRotatableBonds(const Molecule& mol) {
  Cycles cycleData = mol.graph().cycles();

  std::unordered_map<BondIndex, unsigned, boost::hash<BondIndex>> smallestCycle;

  for(const auto& cycleEdges : cycleData) {
    const unsigned cycleSize = cycleEdges.size();

    for(const BondIndex& edge : cycleEdges) {
      auto findIter = smallestCycle.find(edge);

      if(findIter != std::end(smallestCycle)) {
        if(cycleSize < findIter->second) {
          findIter->second = cycleSize;
        }
      } else {
        smallestCycle.emplace(edge, cycleSize);
      }
    }
  }

  double count = 0;
  for(const auto& edge : mol.graph().bonds()) {
    // If the bond is not Single, it cannot be a rotatable bond
    if(mol.graph().bondType(edge) != BondType::Single) {
      continue;
    }

    // If either atom on the bond is terminal, it cannot be rotatable
    if(
      mol.graph().degree(edge.first) == 1
      || mol.graph().degree(edge.second) == 1
    ) {
      continue;
    }

    // If there is an assigned stereopermutator on the edge, it cannot be rotatable
    auto bondStereopermutatorOption = mol.stereopermutators().option(edge);
    if(
      bondStereopermutatorOption
      && bondStereopermutatorOption->assigned() != boost::none
    ) {
      continue;
    }

    auto findIter = smallestCycle.find(edge);

    // Is the edge part of a cycle?
    if(findIter == std::end(smallestCycle)) {
      // If not, the edge counts as a whole rotatable bond
      count += 1.0;
    } else {
      /* Otherwise, the contribution to the number of rotatable bonds is
       * calculated as
       *
       *   max(0.0, (cycle size - 3) / cycle size)
       *
       */
      unsigned cycleSize = findIter->second;
      count += std::max(0.0, (cycleSize - 3.0) / cycleSize);
    }
  }

  return static_cast<unsigned>(
    std::round(count)
  );
}

std::vector<unsigned> rankingEquivalentGroups(const Molecule& mol) {
  PrivateGraph relationshipGraph(mol.graph().V());

  for(
    const AtomStereopermutator& permutator:
    mol.stereopermutators().atomStereopermutators()
  ) {
    const AtomIndex placement = permutator.placement();

    // Skip equivalently ranked substituents of stereogenic permutators
    if(permutator.numAssignments() > 1) {
      continue;
    }

    for(
      const auto& equallyRankedSubstituentGroup:
      permutator.getRanking().substituentRanking
    ) {
      // Skip single-vertex groups, no point to processing those
      if(equallyRankedSubstituentGroup.size() == 1) {
        continue;
      }

      const auto descendants = GraphAlgorithms::bfsUniqueDescendants(
        placement,
        equallyRankedSubstituentGroup,
        mol.graph().inner()
      );

      const auto subgraphs = Temple::map(
        equallyRankedSubstituentGroup,
        [&](const AtomIndex i) -> DescendantSubgraph {
          return DescendantSubgraph(i, descendants, mol.graph().inner());
        }
      );

      // Sequential pairs is enough, no need for all pairs
      Temple::forEach(
        Temple::Adaptors::sequentialPairs(subgraphs),
        [&](const DescendantSubgraph& lhs, const DescendantSubgraph& rhs) {
          /* Avoid an isomorphism if the subgraph is a single vertex, not least
           * because boost::isomorphism yields out of range values in the index
           * map.
           */
          if(lhs.subgraph.V() == 1) {
            const AtomIndex source = lhs.indexMap.left.begin()->first;
            const AtomIndex target = rhs.indexMap.left.begin()->first;
            if(!relationshipGraph.edgeOption(source, target)) {
              relationshipGraph.addEdge(source, target, BondType::Single);
            }
            return;
          }

          const auto isomorphismOption = lhs.subgraph.modularIsomorphism(rhs.subgraph);
          if(!isomorphismOption) {
            throw std::runtime_error("Found no isomorphism for ranking-identical branch subgraphs!");
          }
          const auto isomorphism = isomorphismOption.value();
          const AtomIndex subgraphN = lhs.subgraph.V();
          for(AtomIndex lhsSubgraphV = 0; lhsSubgraphV < subgraphN; ++lhsSubgraphV) {
            const AtomIndex rhsSubgraphV = isomorphism.at(lhsSubgraphV);
            const AtomIndex source = lhs.indexMap.right.at(lhsSubgraphV);
            const AtomIndex target = rhs.indexMap.right.at(rhsSubgraphV);
            if(!relationshipGraph.edgeOption(source, target)) {
              relationshipGraph.addEdge(source, target, BondType::Single);
            }
          }
        }
      );
    }
  }

  /* Transform the relationship graph into connected components and select the
   * first vertex from each component
   */
  std::vector<unsigned> components;
  relationshipGraph.connectedComponents(components);
  return components;
}

std::vector<AtomIndex> rankingDistinctAtoms(const Molecule& mol) {
  const auto components = rankingEquivalentGroups(mol);
  assert(!components.empty());

  // Find the smallest AtomIndex belonging to each of the connected components
  std::unordered_set<unsigned> componentsFound;
  std::vector<AtomIndex> distinct;
  const auto begin = std::begin(components);
  const auto end = std::end(components);
  for(auto iter = begin; iter != end; ++iter) {
    if(componentsFound.count(*iter) == 0) {
      componentsFound.insert(*iter);
      distinct.push_back(iter - begin);
    }
  }

  return distinct;
}

} // namespace Molassembler
} // namespace Scine
