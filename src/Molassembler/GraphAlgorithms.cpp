/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/GraphAlgorithms.h"

#include "Molassembler/Graph.h"
#include "Molassembler/Graph/GraphAlgorithms.h"

/* GraphAlgorithms.h is the public API point for algorithms that act on
 * Graph. In Graph/GraphAlgorithms.h, those algorithms are implemented on
 * the PrivateGraph, which is not accessible using the public API.
 */

namespace Scine {

namespace Molassembler {

std::vector<unsigned> distance(AtomIndex i, const Graph& graph) {
  if(i > graph.N()) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  return GraphAlgorithms::distance(i, graph.inner());
}

} // namespace Molassembler

} // namespace Scine
