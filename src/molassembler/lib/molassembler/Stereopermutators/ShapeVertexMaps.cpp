/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/Stereopermutators/ShapeVertexMaps.h"

#include "molassembler/Stereopermutators/AbstractPermutations.h"

#include "stereopermutation/Stereopermutation.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"

namespace Scine {
namespace molassembler {

SiteToShapeVertexMap siteToShapeVertexMap(
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
  > availableShapeVertices;

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

    availableShapeVertices.emplace_back(
      std::move(positions)
    );
  }

  // For linked sites, we need to find a shape position to place them
  auto placeAndMark = [&](const unsigned shapeVertex) {
    char priority = stereopermutation.characters.at(shapeVertex);

    const unsigned countUpToPosition = std::count(
      std::begin(stereopermutation.characters),
      std::begin(stereopermutation.characters) + shapeVertex,
      priority
    );

    const unsigned correspondingSite = canonicalSites.at(priority - 'A').at(countUpToPosition);

    if(positionMap.at(correspondingSite) == placeholder) {
      unsigned newShapeVertex = availableShapeVertices.at(priority - 'A').front();

      positionMap.at(correspondingSite) = newShapeVertex;

      availableShapeVertices.at(priority - 'A').erase(
        availableShapeVertices.at(priority - 'A').begin()
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

      const unsigned shapeVertex = availableShapeVertices.at(priorityChar - 'A').front();

      availableShapeVertices.at(priorityChar - 'A').erase(
        std::begin(availableShapeVertices.at(priorityChar - 'A'))
      );

      positionMap.at(correspondingSite) = shapeVertex;

      *unusedIndexIter = true;
    }
  }

  // Ensure no shape positions are marked with placeholders
  assert(
    temple::all_of(
      positionMap,
      [](const unsigned shapeVertex) -> bool {
        return shapeVertex != placeholder;
      }
    ) && "A shape position is still marked with a placeholder!"
  );

  return SiteToShapeVertexMap(positionMap);
}

temple::StrongIndexFlatMap<shapes::Vertex, SiteIndex> shapeVertexToSiteIndexMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
) {
  const auto base = siteToShapeVertexMap(stereopermutation, canonicalSites);
  return base.invert();
}

stereopermutation::Stereopermutation stereopermutationFromSiteToShapeVertexMap(
  const SiteToShapeVertexMap& siteToShapeVertexMap,
  const std::vector<RankingInformation::Link>& links,
  const RankingInformation::RankedSitesType& canonicalSites
) {
  /* We need to transform the sites to their canonical ranking characters
   * using the canonical sites and then place them at their shape vertices
   * according to the passed map.
   *
   * Then we need to map the links.
   *
   * From that we should be able to construct a Stereopermutation.
   */
  const unsigned S = siteToShapeVertexMap.size();
  std::vector<char> characters(S);
  for(SiteIndex siteIndex {0}; siteIndex < S; ++siteIndex) {
    auto findIter = std::find_if(
      std::begin(canonicalSites),
      std::end(canonicalSites),
      [siteIndex](const auto& equalSet) -> bool {
        return std::find(
          std::begin(equalSet),
          std::end(equalSet),
          siteIndex
        ) != std::end(equalSet);
      }
    );

    assert(findIter != std::end(canonicalSites));
    unsigned canonicalSetIndex = findIter - std::begin(canonicalSites);
    characters.at(
      siteToShapeVertexMap.at(siteIndex)
    ) = ('A' + canonicalSetIndex);
  }

  // Now for the links.
  stereopermutation::Stereopermutation::OrderedLinks selfReferentialLinks;
  for(auto linkInformation: links) {
    selfReferentialLinks.emplace_back(
      siteToShapeVertexMap.at(linkInformation.sites.first),
      siteToShapeVertexMap.at(linkInformation.sites.second)
    );
  }

  return stereopermutation::Stereopermutation {characters, selfReferentialLinks};
}

} // namespace molassembler
} // namespace Scine
