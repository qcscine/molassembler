/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_CANONICALIZATION_H
#define INCLUDE_MOLASSEMBLER_GRAPH_CANONICALIZATION_H

#include "boost/multiprecision/cpp_int.hpp"
#include "molassembler/Types.h"
#include <vector>

namespace Scine {
namespace molassembler {

class InnerGraph;

namespace hashes {
using WideHashType = boost::multiprecision::uint128_t;
} // namespace hashes

/**
 * @brief Calculate the canonical labeling of a molecule from a coloring
 *   specified by a set of hashes
 *
 * @param inner The inner graph representation of a Molecule
 * @param hashes A flat map of hashes for each vertex
 *
 * @throws std::domain_error If the size of mol's graph exceeds the maximum
 *   value of int. This is the limit for the underlying labeling algorithm.
 *
 * @throws std::invalid_argument If vector of hashes length does not match
 *   the molecule's number of vertices.
 *
 * @return The canonical labeling sequence of atom indices into @p mol.
 */
std::vector<int> canonicalAutomorphism(
  const InnerGraph& inner,
  const std::vector<hashes::WideHashType>& hashes
);

/**
 * @brief Compares two graphs for isomorphism.
 *
 * @param firstGraph The graph of the first molecule to compare
 * @param firstHashes Vertex coloring hashes for the first molecule
 * @param secondGraph The graph of the second molecule to compare
 * @param secondHashes Vertex coloring hashes for the second molecule
 *
 * @return Whether the two graphs are isomorphic
 */
bool isomorphic(
  const InnerGraph& firstGraph,
  const std::vector<hashes::WideHashType>& firstHashes,
  const InnerGraph& secondGraph,
  const std::vector<hashes::WideHashType>& secondHashes
);

} // namespace molassembler
} // namespace Scine

#endif
