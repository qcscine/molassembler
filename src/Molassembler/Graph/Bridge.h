/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Graph and PrivateGraph vertex- & edge descriptor conversions
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H
#define INCLUDE_MOLASSEMBLER_GRAPH_BRIDGE_H

#include "Molassembler/Graph.h"
#include "Molassembler/Graph/PrivateGraph.h"

namespace Scine {
namespace Molassembler {

//! Transform BondIndex to PrivateGraph::Edge
inline PrivateGraph::Edge toInner(const BondIndex& bondIndex, const PrivateGraph& graph) {
  return graph.edge(bondIndex.first, bondIndex.second);
}

//! Transform PrivateGraph::Edge to BondIndex
inline BondIndex toOuter(const PrivateGraph::Edge& edge, const PrivateGraph& graph) {
  return { graph.source(edge), graph.target(edge) };
}

} // namespace Molassembler
} // namespace Scine

#endif
