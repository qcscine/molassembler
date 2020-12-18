/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/GraphAlgorithms.h"

#include "Molassembler/Graph.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Graph/EditDistance.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Subgraphs.h"

// TODO tmp
#include <fstream>
#include <iostream>

/* GraphAlgorithms.h is the public API point for algorithms that act on
 * Graph. In Graph/GraphAlgorithms.h, those algorithms are implemented on
 * the PrivateGraph, which is not accessible using the public API.
 */

namespace Scine {
namespace Molassembler {
namespace {

using UnorderedIndexMap = std::unordered_map<AtomIndex, MultiEdits::ComponentIndexPair>;
std::pair<PrivateGraph, UnorderedIndexMap> condense(const GraphList& list) {
  std::pair<PrivateGraph, UnorderedIndexMap> p;
  unsigned component = 0;
  for(const auto& graphRef : list) {
    const auto indexMap = p.first.merge(graphRef.get().inner());
    for(auto pair : indexMap) {
      p.second.emplace(pair.second, std::make_pair(component, pair.first));
    }
    ++component;
  }
  assert(p.second.size() == p.first.V());
  return p;
}

Edits minimalEdits(
  const PrivateGraph& a,
  const PrivateGraph& b,
  std::unique_ptr<EditCost> cost,
  const Subgraphs::IndexMap& preconditioning = {}
) {
  GraphAlgorithms::EditDistanceForest forest {a, b, std::move(cost), preconditioning};

  Edits edits;
  edits.distance = forest.g[forest.result].costSum;
  const auto traversal = forest.traverse(forest.result);
  edits.indexMap = traversal.bVertices;
  std::reverse(std::begin(edits.indexMap), std::end(edits.indexMap));
  return edits;
}

} // namespace

std::vector<unsigned> distance(AtomIndex source, const Graph& graph) {
  if(source > graph.V()) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  return GraphAlgorithms::distance(source, graph.inner());
}

PredecessorMap shortestPaths(AtomIndex source, const Graph& graph) {
  if(source > graph.V()) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  return PredecessorMap {
    GraphAlgorithms::shortestPaths(source, graph.inner())
  };
}

std::vector<AtomIndex> PredecessorMap::path(const AtomIndex target) const {
  std::vector<AtomIndex> pathVertices;
  AtomIndex position = target;
  while(predecessors.at(position) != position) {
    pathVertices.push_back(position);
    position = predecessors.at(position);
  }
  pathVertices.push_back(position);

  std::reverse(
    std::begin(pathVertices),
    std::end(pathVertices)
  );

  return pathVertices;
}

Edits minimalEdits(const Graph& a, const Graph& b, std::unique_ptr<EditCost> cost) {
  return minimalEdits(a.inner(), b.inner(), std::move(cost));
}

MultiEdits reactionEdits(const GraphList& lhsGraphs, const GraphList& rhsGraphs) {
  const auto lhs = condense(lhsGraphs);
  const auto rhs = condense(rhsGraphs);

  std::ofstream lhsfile("lhs.dot");
  lhsfile << lhs.first.graphviz();
  std::ofstream rhsfile("rhs.dot");
  rhsfile << rhs.first.graphviz();

  // Check preconditions
  auto lhsElements = lhs.first.elementCollection();
  auto rhsElements = rhs.first.elementCollection();
  if(lhsElements.size() != rhsElements.size()) {
    throw std::logic_error("Unequal atom count in reaction. Playing at alchemy?");
  }
  Temple::sort(lhsElements);
  Temple::sort(rhsElements);
  if(lhsElements != rhsElements) {
    throw std::logic_error("Element composition of reaction sides different. Playing at alchemy?");
  }

  /* Precondition the graph edit distance algorithm with maximum common
   * subgraph matches
   */
  const auto edits = Temple::accumulate(
    Subgraphs::maximum(lhs.first, rhs.first),
    Edits { std::numeric_limits<unsigned>::max(), {}},
    [&](Edits carry, const auto& preconditioning) -> Edits {
      std::cout << "MCS precondition size: " << preconditioning.size() << "\n";
      auto preconditionedEdits = minimalEdits(
        lhs.first,
        rhs.first,
        std::make_unique<ElementsConservedCost>(),
        preconditioning
      );

      if(preconditionedEdits.distance < carry.distance) {
        return preconditionedEdits;
      }

      return carry;
    }
  );

  MultiEdits multi;
  multi.distance = edits.distance;

  const unsigned V = edits.indexMap.size();
  assert(V == lhs.first.V());
  for(unsigned i = 0; i < V; ++i) {
    // NOTE: This just assumes there are no vertex epsilons due to the high cost
    multi.indexMap.emplace(
      lhs.second.at(i),
      rhs.second.at(edits.indexMap.at(i))
    );
  }
  return multi;
}

} // namespace Molassembler
} // namespace Scine
