/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/Stereopermutators/SymmetryPositionMaps.h"

#include "stereopermutation/Stereopermutation.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"

namespace Scine {
namespace molassembler {

std::vector<unsigned> siteToSymmetryPositionMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
) {
  /* We are given an stereopermutation within some shape:
   *   Characters: ABADCB
   *   Links: (0, 1), (2, 5)
   * and canonical sites (ranked sets of site indices re-sorted by
   * decreasing set size): {{0, 4}, {2, 1}, {5}, {3}}
   *
   * We need to distribute the site indices (from the canonical sites) to
   * the shape positions (defined by the stereopermutation) that they match (via
   * their ranking character).
   *
   * Additionally, we need to ensure that:
   * - AAAAAA {0, 1}, {2, 3}, {4, 5}
   * - AAAAAA {0, 1}, {2, 4}, {3, 5}
   * have different shape position maps.
   *
   * Link indices (from assignments) specify shape positions that are linked.
   * Since shape positions are NOT exchangeable as two site indices are
   * that rank equally, we need to distribute linked shape positions first,
   * and then distribute the remaining characters afterwards.
   */

  constexpr unsigned placeholder = std::numeric_limits<unsigned>::max();
  const unsigned S = stereopermutation.characters.size();

  std::vector<unsigned> positionMap (S, placeholder);

  /* Process the links */

  /* For every site index within the group of sites of equal priority, we
   * have to keep information on which have been used and which haven't.
   */
  auto usedLists = temple::map(
    canonicalSites,
    [](const auto& equalPrioritySet) -> std::vector<bool> {
      return std::vector<bool>(equalPrioritySet.size(), false);
    }
  );

  /* Additionally, for each canonical character, a limited set of shape
   * positions are available: Those where the passed stereopermutation's characters
   * match the character.
   */
  std::vector<
    std::vector<unsigned>
  > availableSymmetryPositions;

  const char maxChar = temple::max(stereopermutation.characters);
  // For each ranking character
  for(char i = 'A'; i <= maxChar; ++i) {
    std::vector<unsigned> positions;

    // Go through the shape positions
    for(unsigned s = 0; s < S; ++s) {
      if(stereopermutation.characters.at(s) == i) {
        positions.push_back(s);
      }
    }

    availableSymmetryPositions.emplace_back(
      std::move(positions)
    );
  }

  // For linked sites, we need to find a shape position to place them
  auto placeAndMark = [&](const unsigned symmetryPosition) {
    char priority = stereopermutation.characters.at(symmetryPosition);

    unsigned countUpToPosition = std::count(
      std::begin(stereopermutation.characters),
      std::begin(stereopermutation.characters) + symmetryPosition,
      priority
    );

    unsigned correspondingSite = canonicalSites.at(priority - 'A').at(countUpToPosition);

    if(positionMap.at(correspondingSite) == placeholder) {
      unsigned newSymmetryPosition = availableSymmetryPositions.at(priority - 'A').front();

      positionMap.at(correspondingSite) = newSymmetryPosition;

      availableSymmetryPositions.at(priority - 'A').erase(
        availableSymmetryPositions.at(priority - 'A').begin()
      );

      usedLists.at(priority - 'A').at(countUpToPosition) = true;
    }
  };

  // Place all linked indices
  for(const auto& link: stereopermutation.links) {
    placeAndMark(link.first);
    placeAndMark(link.second);
  }

  // Next, process all characters
  for(const auto& priorityChar : stereopermutation.characters) {
    // Get an unused site index for that priority
    const auto unusedIndexIter = std::find(
      std::begin(usedLists.at(priorityChar - 'A')),
      std::end(usedLists.at(priorityChar - 'A')),
      false
    );

    if(unusedIndexIter != std::end(usedLists.at(priorityChar - 'A'))) {
      const unsigned correspondingSite = canonicalSites.at(priorityChar - 'A').at(
        unusedIndexIter - std::begin(usedLists.at(priorityChar - 'A'))
      );

      assert(positionMap.at(correspondingSite) == placeholder);

      const unsigned symmetryPosition = availableSymmetryPositions.at(priorityChar - 'A').front();

      availableSymmetryPositions.at(priorityChar - 'A').erase(
        std::begin(availableSymmetryPositions.at(priorityChar - 'A'))
      );

      positionMap.at(correspondingSite) = symmetryPosition;

      *unusedIndexIter = true;
    }
  }

  // Ensure no shape positions are marked with placeholders
  assert(
    temple::all_of(
      positionMap,
      [](const unsigned symmetryPosition) -> bool {
        return symmetryPosition != placeholder;
      }
    ) && "A shape position is still marked with a placeholder!"
  );

  return positionMap;
}

std::vector<unsigned> symmetryPositionToSiteMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
) {
  auto base = siteToSymmetryPositionMap(stereopermutation, canonicalSites);

  std::vector<unsigned> inverseMap (base.size());

  for(unsigned i = 0; i < base.size(); ++i) {
    inverseMap.at(base.at(i)) = i;
  }

  return inverseMap;
}

} // namespace molassembler
} // namespace Scine
