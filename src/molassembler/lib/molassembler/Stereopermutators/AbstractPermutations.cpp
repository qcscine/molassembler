/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/Stereopermutators/AbstractPermutations.h"

#include <algorithm>

namespace Scine {

namespace molassembler {

RankingInformation::RankedSitesType AbstractStereopermutations::canonicalize(
  RankingInformation::RankedSitesType rankedSites
) {
  std::stable_sort(
    rankedSites.begin(),
    rankedSites.end(),
    [](const auto& setA, const auto& setB) -> bool {
      // Inverted comparison so that larger sets come first
      return setA.size() > setB.size();
    }
  );

  return rankedSites;
}

// Transform canonical ranked sites to canonical characters
std::vector<char> AbstractStereopermutations::transferToSymbolicCharacters(
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

stereopermutation::Stereopermutation::LinksSetType
AbstractStereopermutations::selfReferentialTransform(
  const std::vector<LinkInformation>& rankingLinks,
  const RankingInformation::RankedSitesType& canonicalSites
) {
  stereopermutation::Stereopermutation::LinksSetType links;

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

    links.emplace(
      std::min(a, b),
      std::max(a, b)
    );
  }

  return links;
}

std::vector<char> AbstractStereopermutations::makeStereopermutationCharacters(
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

  for(const auto& index : sitesAtShapeVertices) {
    const auto findIter = std::find(
      flattenedIndices.begin(),
      flattenedIndices.end(),
      index
    );

    newStereopermutationCharacters.push_back(
      canonicalStereopermutationCharacters.at(
        findIter - flattenedIndices.begin()
      )
    );
  }

  return newStereopermutationCharacters;
}

AbstractStereopermutations::AbstractStereopermutations(
  const RankingInformation& ranking,
  const Symmetry::Shape shape
) : canonicalSites(canonicalize(ranking.siteRanking)),
    symbolicCharacters(transferToSymbolicCharacters(canonicalSites)),
    selfReferentialLinks(selfReferentialTransform(ranking.links, canonicalSites)),
    permutations(
      stereopermutation::uniquesWithWeights(
        stereopermutation::Stereopermutation {
          symbolicCharacters,
          selfReferentialLinks
        },
        shape,
        false
      )
    )
{}

} // namespace molassembler

} // namespace Scine
