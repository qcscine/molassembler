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
Stereopermutations::Stereopermutation::Occupation
Abstract::transferToOccupation(
  const RankingInformation::RankedSitesType& canonicalSites
) {
  std::vector<unsigned> rank;
  unsigned currentRank = 0;
  for(const auto& equalPrioritySet : canonicalSites) {
    for(unsigned i = 0; i < equalPrioritySet.size(); ++i) {
      rank.push_back(currentRank);
    }

    ++currentRank;
  }

  return Stereopermutations::Stereopermutation::Occupation {std::move(rank)};
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

  const auto getRankedVertex = [&canonicalSites](const SiteIndex site) -> Shapes::Vertex {
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

  using Link = Stereopermutations::Stereopermutation::Link;

  return Temple::sorted(
    Temple::map(
      rankingLinks,
      [&](const auto& siteLink) -> Link {
        Link vertexLink {
          getRankedVertex(siteLink.sites.first),
          getRankedVertex(siteLink.sites.second)
        };

        if(vertexLink.first > vertexLink.second) {
          std::swap(vertexLink.first, vertexLink.second);
        }

        if(vertexLink.first == vertexLink.second) {
          throw std::runtime_error("Got same ranked vertex");
        }

        return vertexLink;
      }
    )
  );
}

Stereopermutations::Stereopermutation::Occupation
Abstract::makeOccupation(
  const RankingInformation::RankedSitesType& canonicalSites,
  const Stereopermutations::Stereopermutation::Occupation& canonicalOccupation,
  const Temple::StrongIndexPermutation<Shapes::Vertex, SiteIndex>& sitesAtShapeVertices
) {
  // Replace the site indices by their new ranking characters
  std::vector<SiteIndex> flattenedSiteIndices;
  for(const auto& equalPrioritySet : canonicalSites) {
    for(const auto& index : equalPrioritySet) {
      flattenedSiteIndices.push_back(index);
    }
  }

  std::vector<unsigned> occupation;

  for(const auto keyValuePair : sitesAtShapeVertices) {
    const auto findIter = std::find(
      std::begin(flattenedSiteIndices),
      std::end(flattenedSiteIndices),
      keyValuePair.second
    );

    assert(findIter != std::end(flattenedSiteIndices));

    occupation.push_back(
      canonicalOccupation.at(
        Shapes::Vertex {
          static_cast<unsigned>(findIter - flattenedSiteIndices.begin())
        }
      )
    );
  }

  return Stereopermutations::Stereopermutation::Occupation {occupation};
}

Abstract::Abstract(
  const RankingInformation& ranking,
  const Shapes::Shape shape
) : canonicalSites(canonicalize(ranking.siteRanking)),
    occupation(transferToOccupation(canonicalSites)),
    selfReferentialLinks(selfReferentialTransform(ranking.links, canonicalSites)),
    permutations(
      Stereopermutations::uniques(
        Stereopermutations::Stereopermutation {occupation, selfReferentialLinks},
        shape,
        false
      )
    )
{}

} // namespace Stereopermutators
} // namespace Molassembler
} // namespace Scine
