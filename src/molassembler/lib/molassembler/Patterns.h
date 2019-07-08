/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Molecule patterns for better rapid prototyping
 */

#ifndef INCLUDE_MOLASSEMBLER_PATTERNS_H
#define INCLUDE_MOLASSEMBLER_PATTERNS_H

#include "molassembler/Types.h"
#include <vector>

namespace Scine {
namespace molassembler {

// Forward-declarations
class Molecule;

namespace patterns {

using PlugType = std::pair<Molecule, std::vector<AtomIndex>>;

PlugType methyl();

} // namespace patterns
} // namespace molassembler
} // namespace Scine

#endif
