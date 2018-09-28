// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/GraphAlgorithms.h"

#include "molassembler/OuterGraph.h"
#include "molassembler/Graph/GraphAlgorithms.h"

/* GraphAlgorithms.h is the public API point for algorithms that act on
 * OuterGraph. In Graph/GraphAlgorithms.h, those algorithms are implemented on
 * the InnerGraph, which is not accessible using the public API.
 */

namespace molassembler {

std::vector<unsigned> distance(AtomIndex i, const OuterGraph& graph) {
  if(i > graph.N()) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  return GraphAlgorithms::distance(i, graph.inner());
}

} // namespace molassembler
