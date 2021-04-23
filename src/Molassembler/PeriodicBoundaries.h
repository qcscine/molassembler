/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Extras for dealing with periodic boundaries
 */

#ifndef INCLUDE_MOLASSEMBLER_PERIODIC_BOUNDARIES_H
#define INCLUDE_MOLASSEMBLER_PERIODIC_BOUNDARIES_H

#include "Molassembler/Types.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace Scine {
namespace Molassembler {

class Graph;

struct PeriodicBoundaryDuplicates {
  std::unordered_set<AtomIndex> uninterestingAtoms;
  std::unordered_map<AtomIndex, AtomIndex> ghostAtomMap;
};

struct SubstitutionsGenerator {
  //! List of real-space atom index to ghost atom index replacements
  using SubstitutionList = std::vector<std::pair<AtomIndex, AtomIndex>>;

  //! Map from interesting atom to list of real-to-ghost atom replacements
  using SubstitutionMap = std::unordered_map<AtomIndex, SubstitutionList>;

  /*! @brief Remove edges from interesting real-space atoms to ghost-space atoms
   */
  static SubstitutionMap removeGhosts(
    Graph& graph,
    const PeriodicBoundaryDuplicates& periodics
  );
};

} // namespace Molassembler
} // namespace Scine

#endif
