/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "Molassembler/Stereopermutators/AbstractPermutations.h"

#include "Molassembler/Temple/Functional.h"
#include "boost/optional.hpp"

#include <algorithm>
#include <cassert>

namespace Scine {
namespace Molassembler {
namespace Stereopermutators {

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

Stereopermutations::Stereopermutation::OrderedLinks
Abstract::selfReferentialTransform(
  const std::vector<RankingInformation::Link>& rankingLinks,
  const RankingInformation::RankedSitesType& canonicalSites
) {
  if(
    Temple::any_of(
      rankingLinks,
      [](auto&& l) { return l.sites.first == l.sites.second; }
    )
  ) {
    throw std::out_of_range("Links are invalid");
  }

  auto getRankedVertex = [&canonicalSites](const SiteIndex site) -> Shapes::Vertex {
    Shapes::Vertex vertex(0);
    for(const auto& equalSitesSet : canonicalSites) {
      for(const SiteIndex& rankedSiteIndex : equalSitesSet) {
        if(rankedSiteIndex == site) {
          return vertex;
        }

        ++vertex;
      }
    }

    throw std::logic_error("Site index not found in ranked sites");
  };

  return Temple::sorted(
    Temple::map(
      rankingLinks,
      [&](const auto& link) -> Stereopermutations::Stereopermutation::Link {
        auto firstVertex = getRankedVertex(link.sites.first);
        auto secondVertex = getRankedVertex(link.sites.second);

        if(firstVertex > secondVertex) {
          std::swap(firstVertex, secondVertex);
        }

        if(firstVertex == secondVertex) {
          throw std::runtime_error("Got same ranked vertex");
        }

        return {firstVertex, secondVertex};
      }
    )
  );
}

std::vector<char> Abstract::makeStereopermutationCharacters(
  const RankingInformation::RankedSitesType& canonicalSites,
  const std::vector<char>& canonicalStereopermutationCharacters,
  const Temple::StrongIndexFlatMap<Shapes::Vertex, SiteIndex>& sitesAtShapeVertices
) {
  // Replace the site indices by their new ranking characters
  std::vector<SiteIndex> flattenedIndices;
  for(const auto& equalPrioritySet : canonicalSites) {
    for(const auto& index : equalPrioritySet) {
      flattenedIndices.push_back(index);
    }
  }

  std::vector<char> newStereopermutationCharacters;

  for(const SiteIndex siteIndex : sitesAtShapeVertices) {
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
  const Shapes::Shape shape
) : canonicalSites(canonicalize(ranking.siteRanking)),
    symbolicCharacters(transferToSymbolicCharacters(canonicalSites)),
    selfReferentialLinks(selfReferentialTransform(ranking.links, canonicalSites)),
    permutations(
      Stereopermutations::uniques(
        Stereopermutations::Stereopermutation {
          symbolicCharacters,
          selfReferentialLinks
        },
        shape,
        false
      )
    )
{}

} // namespace Stereopermutators
} // namespace Molassembler
} // namespace Scine
