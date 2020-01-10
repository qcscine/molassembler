/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/Stereopermutators/AbstractPermutations.h"

#include <algorithm>
#include <cassert>

namespace Scine {
namespace molassembler {
namespace stereopermutators {

RankingInformation::RankedSitesType Abstract::canonicalize(
  RankingInformation::RankedSitesType rankedSites
) {
  std::stable_sort(
    std::begin(rankedSites),
    std::end(rankedSites),
    [](const auto& setA, const auto& setB) -> bool {
      // Inverted comparison so that larger sets come first
      return setA.size() > setB.size();
    }
  );

  return rankedSites;
}

// Transform canonical ranked sites to canonical characters
std::vector<char> Abstract::transferToSymbolicCharacters(
  const RankingInformation::RankedSitesType& canonicalSites
) {
  std::vector<char> characters;

  char currentChar = 'A';
  for(const auto& equalPrioritySet : canonicalSites) {
    for(unsigned i = 0; i < equalPrioritySet.size(); ++i) {
      characters.push_back(currentChar);
    }

    ++currentChar;
  }

  return characters;
}

stereopermutation::Stereopermutation::OrderedLinks
Abstract::selfReferentialTransform(
  const std::vector<LinkInformation>& rankingLinks,
  const RankingInformation::RankedSitesType& canonicalSites
) {
  stereopermutation::Stereopermutation::OrderedLinks links;

  for(const auto& link : rankingLinks) {
    auto getRankedPosition = [&canonicalSites](const unsigned siteIndex) -> unsigned {
      unsigned position = 0;
      for(const auto& equalSitesSet : canonicalSites) {
        for(const auto& rankedSiteIndex : equalSitesSet) {
          if(rankedSiteIndex == siteIndex) {
            return position;
          }

          ++position;
        }
      }

      throw std::logic_error("Site index not found in ranked sites");
    };

    const unsigned a = getRankedPosition(link.indexPair.first);
    const unsigned b = getRankedPosition(link.indexPair.second);

    links.emplace_back(
      std::min(a, b),
      std::max(a, b)
    );
  }

  std::sort(std::begin(links), std::end(links));
  return links;
}

std::vector<char> Abstract::makeStereopermutationCharacters(
  const RankingInformation::RankedSitesType& canonicalSites,
  const std::vector<char>& canonicalStereopermutationCharacters,
  const std::vector<unsigned>& sitesAtShapeVertices
) {
  // Replace the site indices by their new ranking characters
  std::vector<unsigned> flattenedIndices;
  for(const auto& equalPrioritySet : canonicalSites) {
    for(const auto& index : equalPrioritySet) {
      flattenedIndices.push_back(index);
    }
  }

  std::vector<char> newStereopermutationCharacters;

  for(const auto& siteIndex : sitesAtShapeVertices) {
    const auto findIter = std::find(
      std::begin(flattenedIndices),
      std::end(flattenedIndices),
      siteIndex
    );

    assert(findIter != std::end(flattenedIndices));

    newStereopermutationCharacters.push_back(
      canonicalStereopermutationCharacters.at(
        findIter - flattenedIndices.begin()
      )
    );
  }

  return newStereopermutationCharacters;
}

Abstract::Abstract(
  const RankingInformation& ranking,
  const shapes::Shape shape
) : canonicalSites(canonicalize(ranking.siteRanking)),
    symbolicCharacters(transferToSymbolicCharacters(canonicalSites)),
    selfReferentialLinks(selfReferentialTransform(ranking.links, canonicalSites)),
    permutations(
      stereopermutation::uniques(
        stereopermutation::Stereopermutation {
          symbolicCharacters,
          selfReferentialLinks
        },
        shape,
        false
      )
    )
{}

} // namespace stereopermutators
} // namespace molassembler
} // namespace Scine
