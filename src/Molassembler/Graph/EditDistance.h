/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_EDIT_DISTANCE_H
#define INCLUDE_MOLASSEMBLER_GRAPH_EDIT_DISTANCE_H

#include "Molassembler/Graph/PrivateGraph.h"

namespace Scine {
namespace Molassembler {
namespace GraphAlgorithms {

/**
 * @brief Exact graph edit distance calculation
 *
 * Graph edit distance is symmetric, so order of arguments is irrelevant.
 *
 * @param a First graph to calculate edit distance for
 * @param b Second graph to calculate edit distance for
 *
 * @return Graph edit distance metric
 */
unsigned editDistance(const PrivateGraph& a, const PrivateGraph& b);

} // namespace GraphAlgorithms
} // namespace Molassembler
} // namespace Scine

#endif
