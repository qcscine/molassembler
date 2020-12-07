/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/GraphAlgorithms.h"

#include "Molassembler/Graph.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Graph/EditDistance.h"

/* GraphAlgorithms.h is the public API point for algorithms that act on
 * Graph. In Graph/GraphAlgorithms.h, those algorithms are implemented on
 * the PrivateGraph, which is not accessible using the public API.
 */

namespace Scine {
namespace Molassembler {

std::vector<unsigned> distance(AtomIndex source, const Graph& graph) {
  if(source > graph.N()) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  return GraphAlgorithms::distance(source, graph.inner());
}

PredecessorMap shortestPaths(AtomIndex source, const Graph& graph) {
  if(source > graph.N()) {
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

unsigned editDistance(const Graph& a, const Graph& b) {
  return GraphAlgorithms::editDistance(a.inner(), b.inner());
}

} // namespace Molassembler
} // namespace Scine
