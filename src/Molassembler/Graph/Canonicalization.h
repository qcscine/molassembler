/*! @file
 * @brief Provides canonicalization of Molecule instances
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 *
 *
 * Graph canonicalization is an algorithm that generates an isomorphic graph
 * to its input. The important property is that the algorithm's output graph is
 * identical for all graphs that are isomorphic to its input. For two canonical
 * graphs, an isomorphism test reduces to a fast identity test.
 *
 * Naturally, graph canonicalization is computationally at least as hard as an
 * isomorphism test on non-canonical graphs, so little time (if any) is saved
 * by canonicalizing two molecules separately, then using the identity
 * comparison. However, it is worth canonicalizing molecules when it is to be
 * expected that repeated isomorphism computations need to be performed upon
 * many molecules. For instance, checking whether a molecule is already present
 * in a set of molecules can be sped up considerably.
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_CANONICALIZATION_H
#define INCLUDE_MOLASSEMBLER_GRAPH_CANONICALIZATION_H

#include "boost/multiprecision/cpp_int.hpp"
#include "Molassembler/Types.h"
#include <vector>

namespace Scine {
namespace Molassembler {

class PrivateGraph;

namespace Hashes {
using WideHashType = boost::multiprecision::uint128_t;
} // namespace Hashes

/** @brief Calculate the canonical labeling of a molecule from a coloring
 *   specified by a set of hashes
 *
 * @complexity{Underlying algorithm is exponential-time, but we use graph
 * invariants in addition, which should give us sub-exponential complexity for
 * common cases.}
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
 * @warning The canonical form of the molecule is influenced by the provided
 *   hashes. Any comparisons on the canonical form must use the same comparison
 *   bitmask as when generating the hashes for this function call here.
 *
 * @return The canonical labeling sequence of atom indices into @p mol.
 */
std::vector<int> canonicalAutomorphism(
  const PrivateGraph& inner,
  const std::vector<Hashes::WideHashType>& hashes
);

} // namespace Molassembler
} // namespace Scine

#endif
