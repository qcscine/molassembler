/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "Molassembler/Stereopermutators/ShapeVertexMaps.h"

#include "Molassembler/Stereopermutators/AbstractPermutations.h"

#include "Molassembler/Stereopermutation/Stereopermutation.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/constexpr/Numeric.h"
#include "Molassembler/Temple/Poset.h"
#include "Molassembler/Shapes/Properties.h"

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/isomorphism.hpp"

namespace Scine {
namespace Molassembler {
namespace {

using PlainGraphType = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
using Vertex = typename PlainGraphType::vertex_descriptor;

struct SiteIndexColor {
  std::vector<unsigned> rankingColor;

  using argument_type = Vertex;
  using result_type = unsigned;

  inline SiteIndexColor(std::vector<unsigned> colors)
    : rankingColor(std::move(colors)) {}

  inline unsigned operator() (const Vertex i) const {
    return rankingColor.at(i);
  }
};

struct VertexColor {
  const Stereopermutations::Stereopermutation* const stereopermutationPtr;

  using argument_type = Vertex;
  using result_type = unsigned;

  inline VertexColor(const Stereopermutations::Stereopermutation& s)
    : stereopermutationPtr(&s) {}

  inline unsigned operator() (const Vertex i) const {
    return stereopermutationPtr->occupation.at(Shapes::Vertex {static_cast<unsigned>(i)});
  }
};

BOOST_CONCEPT_ASSERT((
  boost::AdaptableUnaryFunctionConcept<SiteIndexColor, unsigned, Vertex>
));
BOOST_CONCEPT_ASSERT((
  boost::AdaptableUnaryFunctionConcept<VertexColor, unsigned, Vertex>
));

template<typename Invariant1, typename Invariant2>
void mapUnmappedVertices(
  Invariant1&& invariant1,
  Invariant2&& invariant2,
  std::vector<Vertex>& mapping
) {
  const unsigned V = mapping.size();
  const Vertex nullVertex = PlainGraphType::null_vertex();
  std::vector<Vertex> unmapped1;
  for(Vertex i = 0; i < V; ++i) {
    if(mapping.at(i) == nullVertex) {
      unmapped1.push_back(i);
    }
  }

  if(unmapped1.empty()) {
    return;
  }

  // Which vertices in g2 are unmapped?
  std::vector<Vertex> mapped2;
  mapped2.reserve(V - unmapped1.size());
  for(Vertex i = 0; i < V; ++i) {
    Vertex mappedVertex = mapping.at(i);
    if(mappedVertex != nullVertex) {
      mapped2.push_back(mappedVertex);
    }
  }
  Temple::sort(mapped2);

  auto allVertices = Temple::iota<Vertex>(V);
  std::vector<Vertex> unmapped2;
  std::set_difference(
    std::begin(allVertices),
    std::end(allVertices),
    std::begin(mapped2),
    std::end(mapped2),
    std::back_inserter(unmapped2)
  );

  using Invariant2Value = typename Invariant2::result_type;
  std::unordered_multimap<Invariant2Value, Vertex> unmatchedMap;

  for(Vertex v : unmapped2) {
    unmatchedMap.emplace(invariant2(v), v);
  }

  // Map the remainder
  for(Vertex v : unmapped1) {
    const auto invariant = invariant1(v);
    const auto findIter = unmatchedMap.find(invariant);
    if(findIter == std::end(unmatchedMap)) {
      throw std::logic_error("Mismatched invariants");
    }
    mapping.at(v) = findIter->second;
    unmatchedMap.erase(findIter);
  }
}

} // namespace

template<typename T>
std::vector<unsigned> builtRankingColors(
  const RankingInformation::RankedSitesType& rankingInformation,
  const std::vector<std::vector<T>>& groups,
  const unsigned maxSites
  ) {
  std::vector<unsigned> grouping(maxSites, 0);
  if(groups.size()) {
    grouping = Temple::map(
      Temple::iota<unsigned>(maxSites),
      [&](unsigned site) -> unsigned {
        const auto findIter = Temple::find_if(
          groups,
          [&](const auto& equalPositions) -> bool {
            return Temple::find(equalPositions, site) != std::end(equalPositions);
          }
        );

        if(findIter == std::end(groups)) {
          throw std::logic_error("Could not find site/vertex in position groups.");
        }

        return findIter - std::begin(groups);
      }
    );
  } // siteGroups.size()

  std::vector<unsigned> rankingColors(maxSites);
  unsigned currentMaxRank = rankingInformation.size() - 1;
  for(unsigned rank = 0; rank < rankingInformation.size(); ++rank) {
    const auto& sitesInRank = rankingInformation[rank];
    if(sitesInRank.size() == 1) {
      const unsigned siteZero = sitesInRank[0];
      rankingColors[siteZero] = rank;
    } else {
      Temple::Poset<SiteIndex> siteIndexPoset {sitesInRank};
      siteIndexPoset.orderUnordered(
        [&](const SiteIndex i, const SiteIndex j) -> bool {
          return grouping[i] < grouping[j];
        }
      );
      const auto& orderedSubset = siteIndexPoset.extract();
      for(const auto& site : orderedSubset[0]) {
        rankingColors[site] = rank;
      } // for site
      for(unsigned i = 1; i < orderedSubset.size(); ++i) {
        ++currentMaxRank;
        for(const auto& site : orderedSubset[i]) {
          rankingColors[site] = currentMaxRank;
        } // for site
      } // for i
    } // else sitesInRank.size() == 1
  }
  return rankingColors;
}

SiteToShapeVertexMap siteToShapeVertexMap(
  const Stereopermutations::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites,
  const std::vector<RankingInformation::Link>& siteLinks,
  std::vector<std::vector<SiteIndex>> siteGroups,
  std::vector<std::vector<Shapes::Vertex>> vertexGroups
) {
  const bool linksExist = (!stereopermutation.links.empty() || !siteLinks.empty());
  const bool mappingIsBijective = (canonicalSites.size() == stereopermutation.occupation.size());

  if(linksExist && !mappingIsBijective) {
    if(stereopermutation.links.size() != siteLinks.size()) {
      throw std::logic_error("Mismatched stereoperm and site link sizes");
    }

    const unsigned S = stereopermutation.occupation.size();
    PlainGraphType siteGraph(S);
    for(const auto& siteLink : siteLinks) {
      boost::add_edge(siteLink.sites.first, siteLink.sites.second, siteGraph);
    }

    PlainGraphType vertexGraph(S);
    for(const auto& vertexLink : stereopermutation.links) {
      boost::add_edge(vertexLink.first, vertexLink.second, vertexGraph);
    }

    if(siteGroups.size() == 0 || vertexGroups.size() == 0) {
      siteGroups.clear();
      vertexGroups.clear();
    }

    /*
     * Build ranking colors for site and vertices.
     * We may include the position groups for the sites and vertices in the
     * color definition for the graph-invariants. In this way, we ensure that
     * any isomorphism of both link-models (vertexGraph/siteGraph) conserve
     * these groups. Otherwise, if we have a degenerate solution, we may get
     * an index to vertex map that violates the groups, and therefore corresponds
     * to a different molecular geometry.
     */
    std::vector<unsigned> siteRankingColors = builtRankingColors<SiteIndex>(
      canonicalSites, siteGroups, S);

    const auto& occupation = stereopermutation.occupation;
    unsigned maxRank = canonicalSites.size();
    RankingInformation::RankedSitesType vertexRanking(maxRank);
    for(unsigned iv = 0; iv < S; ++iv) {
      const auto& rank = occupation.at(Shapes::Vertex {iv});
      vertexRanking.at(rank).emplace_back(iv);
    }

    std::vector<unsigned> vertexRankingColors = builtRankingColors<Shapes::Vertex>(
      vertexRanking, vertexGroups, S);

    std::vector<Vertex> indexMap(S);
    const unsigned maxColor = Temple::max(siteRankingColors) + 1;
    const unsigned maxVertexColor = Temple::max(vertexRankingColors) + 1;
    if(maxColor != maxVertexColor) {
      // Both link-models cannot be isomorphic.
      return SiteToShapeVertexMap{};
    }

    const bool isomorphic = boost::isomorphism(
      siteGraph,
      vertexGraph,
      boost::make_safe_iterator_property_map(
        indexMap.begin(),
        S,
        boost::get(boost::vertex_index, siteGraph)
      ),
      SiteIndexColor(siteRankingColors),
      SiteIndexColor(vertexRankingColors),
      maxColor,
      boost::get(boost::vertex_index, siteGraph),
      boost::get(boost::vertex_index, vertexGraph)
    );

    if(!isomorphic) {
      // Return an empty site to vertex map, if the stereopermutation did not
      // yield a valid isomorphism.
      return SiteToShapeVertexMap{};
    }

    mapUnmappedVertices(
      SiteIndexColor(siteRankingColors),
      SiteIndexColor(vertexRankingColors),
      indexMap
    );

    if(Temple::any_of(indexMap, [S](const Vertex i) { return i >= S; })) {
      throw std::logic_error("Isomorphism index map contains out of bounds vertex indices");
    }

    return SiteToShapeVertexMap::from(
      Temple::map(indexMap, [](auto&& i) { return Shapes::Vertex(i); })
    );
  }

  // If there are no links, this is a little easier
  constexpr Shapes::Vertex placeholder(std::numeric_limits<unsigned>::max());
  const unsigned S = stereopermutation.occupation.size();
  std::vector<Shapes::Vertex> map(S, placeholder);
  /* Maps from site indices to the ranking character (that we can compare with
   * in the stereopermutation characters) by searching for it in the
   * canonicalSites.
   */
  std::vector<Stereopermutations::Rank> siteRanks = Temple::map(
    Temple::iota<SiteIndex>(S),
    [&](SiteIndex site) -> Stereopermutations::Rank {
      const auto findIter = Temple::find_if(
        canonicalSites,
        [&](const auto& equallyRankedSites) -> bool {
          return Temple::find(equallyRankedSites, site) != std::end(equallyRankedSites);
        }
      );

      if(findIter == std::end(canonicalSites)) {
        throw std::logic_error("Could not find site in canonicalSites");
      }

      return Stereopermutations::Rank {
        static_cast<Stereopermutations::Rank::value_type>(findIter - std::begin(canonicalSites))
      };
    }
  );

  auto firstAvailableVertex = [&](const SiteIndex site) -> Shapes::Vertex {
    const auto siteRankedCharacter = siteRanks.at(site);
    for(Shapes::Vertex i {0}; i < S; ++i) {
      if(
        stereopermutation.occupation.at(i) == siteRankedCharacter
        && Temple::find(map, i) == std::end(map)
      ) {
        return Shapes::Vertex {i};
      }
    }

    throw std::logic_error("No available vertex found for this site");
  };

  for(unsigned i = 0; i < S; ++i) {
    if(map.at(i) == placeholder) {
      map.at(i) = firstAvailableVertex(SiteIndex(i));
    }
  }

  return SiteToShapeVertexMap::from(map);
}

Temple::StrongIndexPermutation<Shapes::Vertex, SiteIndex> shapeVertexToSiteIndexMap(
  const Stereopermutations::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites,
  const std::vector<RankingInformation::Link>& siteLinks
) {
  return siteToShapeVertexMap(stereopermutation, canonicalSites, siteLinks).inverse();
}

Stereopermutations::Stereopermutation stereopermutationFromSiteToShapeVertexMap(
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
  std::vector<unsigned> occupation = Temple::map(
    Temple::iota<SiteIndex>(S),
    [&](const SiteIndex siteIndex) -> unsigned {
      return Temple::index_if(
        canonicalSites,
        [siteIndex](const auto& equalSet) -> bool {
          return Temple::makeContainsPredicate(equalSet)(siteIndex);
        }
      );
    }
  );

  occupation = siteToShapeVertexMap.permutation.apply(occupation);

  // Now for the links.
  Stereopermutations::Stereopermutation::OrderedLinks selfReferentialLinks;
  for(const auto& linkInformation: links) {
    selfReferentialLinks.emplace_back(
      siteToShapeVertexMap.at(linkInformation.sites.first),
      siteToShapeVertexMap.at(linkInformation.sites.second)
    );
  }

  return Stereopermutations::Stereopermutation {
    Stereopermutations::Stereopermutation::Occupation {occupation},
    selfReferentialLinks
  };
}

} // namespace Molassembler
} // namespace Scine
