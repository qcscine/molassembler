/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/GraphAlgorithms.h"

#include "molassembler/Graph.h"
#include "molassembler/Graph/GraphAlgorithms.h"

/* GraphAlgorithms.h is the public API point for algorithms that act on
 * Graph. In Graph/GraphAlgorithms.h, those algorithms are implemented on
 * the PrivateGraph, which is not accessible using the public API.
 */

namespace Scine {

namespace molassembler {

std::vector<unsigned> distance(AtomIndex i, const Graph& graph) {
  if(i > graph.N()) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  return graph_algorithms::distance(i, graph.inner());
}

} // namespace molassembler

} // namespace Scine
