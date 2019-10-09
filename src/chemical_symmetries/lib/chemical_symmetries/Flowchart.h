/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Flowchart coordinates for point group symmetries
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_FLOWCHART_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_FLOWCHART_H

#include "chemical_symmetries/Recognition.h"

namespace Scine{
namespace Symmetry {

/**
 * @brief
 *
 * @param normalizedPositions
 *
 * @pre standardizeTop should have been called on normalizedPositions, but not
 * reorientAsymmetricTop (duplicate work otherwise)
 *
 * @return
 */
std::pair<PointGroup, double> flowchart(const PositionCollection& normalizedPositions);

} // namespace Symmetry
} // namespace Scine

#endif
