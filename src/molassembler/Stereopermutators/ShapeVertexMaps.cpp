/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/Stereopermutators/ShapeVertexMaps.h"

#include "molassembler/Stereopermutators/AbstractPermutations.h"

#include "molassembler/Stereopermutation/Stereopermutation.h"
#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/constexpr/Numeric.h"

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/isomorphism.hpp"

namespace Scine {
namespace molassembler {
namespace {

using PlainGraphType = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
using Vertex = typename PlainGraphType::vertex_descriptor;

struct SiteIndexColor {
  std::vector<unsigned> rankingColor;

  using argument_type = Vertex;
  using result_type = unsigned;

  inline SiteIndexColor(std::vector<unsigned> colors)
    : rankingColor(colors) {}

  inline unsigned operator() (const Vertex i) const {
    return rankingColor.at(i);
  }
};

struct VertexColor {
  const stereopermutation::Stereopermutation* const stereopermutationPtr;

  using argument_type = Vertex;
  using result_type = unsigned;

  inline VertexColor(const stereopermutation::Stereopermutation& s)
    : stereopermutationPtr(&s) {}

  inline unsigned operator() (const Vertex i) const {
    return stereopermutationPtr->characters.at(i) - 'A';
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
  temple::inplace::sort(mapped2);

  auto allVertices = temple::iota<Vertex>(V);
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

SiteToShapeVertexMap siteToShapeVertexMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites,
  const std::vector<RankingInformation::Link>& siteLinks
) {
  const bool linksExist = (!stereopermutation.links.empty() || !siteLinks.empty());
  const bool mappingIsBijective = (canonicalSites.size() == stereopermutation.characters.size());

  if(linksExist && !mappingIsBijective) {
    if(stereopermutation.links.size() != siteLinks.size()) {
      throw std::logic_error("Mismatched stereoperm and site link sizes");
    }

    const unsigned S = stereopermutation.characters.size();
    PlainGraphType siteGraph(S);
    for(const auto& siteLink : siteLinks) {
      boost::add_edge(siteLink.sites.first, siteLink.sites.second, siteGraph);
    }

    PlainGraphType vertexGraph(S);
    for(const auto& vertexLink : stereopermutation.links) {
      boost::add_edge(vertexLink.first, vertexLink.second, vertexGraph);
    }

    auto siteRankingColors = temple::map(
      temple::iota<SiteIndex>(S),
      [&](SiteIndex site) -> unsigned {
        const auto findIter = temple::find_if(
          canonicalSites,
          [&](const auto& equallyRankedSites) -> bool {
            return temple::find(equallyRankedSites, site) != std::end(equallyRankedSites);
          }
        );

        if(findIter == std::end(canonicalSites)) {
          throw std::logic_error("Could not find site in canonicalSites");
        }

        return findIter - std::begin(canonicalSites);
      }
    );

    std::vector<Vertex> indexMap(S);
    const unsigned maxColor = temple::max(siteRankingColors) + 1;

    const bool isomorphic = boost::isomorphism(
      siteGraph,
      vertexGraph,
      boost::make_safe_iterator_property_map(
        indexMap.begin(),
        S,
        boost::get(boost::vertex_index, siteGraph)
      ),
      SiteIndexColor(siteRankingColors),
      VertexColor(stereopermutation),
      maxColor,
      boost::get(boost::vertex_index, siteGraph),
      boost::get(boost::vertex_index, vertexGraph)
    );

    if(!isomorphic) {
      throw std::logic_error("Graphs of site and shape vertex links are not isomorphic");
    }

    mapUnmappedVertices(
      SiteIndexColor(siteRankingColors),
      VertexColor(stereopermutation),
      indexMap
    );

    if(temple::any_of(indexMap, [S](const Vertex i) { return i >= S; })) {
      throw std::logic_error("Isomorphism index map contains out of bounds vertex indices");
    }

    return SiteToShapeVertexMap(
      temple::map(indexMap, [](auto&& i) { return shapes::Vertex(i); })
    );
  }

  // If there are no links, this is a little easier
  constexpr shapes::Vertex placeholder(std::numeric_limits<unsigned>::max());
  const unsigned S = stereopermutation.characters.size();
  std::vector<shapes::Vertex> map(S, placeholder);
  /* Maps from site indices to the ranking character (that we can compare with
   * in the stereopermutation characters) by searching for it in the
   * canonicalSites.
   */
  std::vector<char> siteRankingCharacters = temple::map(
    temple::iota<SiteIndex>(S),
    [&](SiteIndex site) -> char {
      const auto findIter = temple::find_if(
        canonicalSites,
        [&](const auto& equallyRankedSites) -> bool {
          return temple::find(equallyRankedSites, site) != std::end(equallyRankedSites);
        }
      );

      if(findIter == std::end(canonicalSites)) {
        throw std::logic_error("Could not find site in canonicalSites");
      }

      return 'A' + (findIter - std::begin(canonicalSites));
    }
  );

  auto firstAvailableVertex = [&](const SiteIndex site) -> shapes::Vertex {
    const char siteRankedCharacter = siteRankingCharacters.at(site);
    for(unsigned i = 0; i < S; ++i) {
      if(
        stereopermutation.characters.at(i) == siteRankedCharacter
        && temple::find(map, shapes::Vertex(i)) == std::end(map)
      ) {
        return shapes::Vertex(i);
      }
    }

    throw std::logic_error("No available vertex found for this site");
  };

  for(unsigned i = 0; i < S; ++i) {
    if(map.at(i) == placeholder) {
      map.at(i) = firstAvailableVertex(SiteIndex(i));
    }
  }

  return SiteToShapeVertexMap(map);
}

temple::StrongIndexFlatMap<shapes::Vertex, SiteIndex> shapeVertexToSiteIndexMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites,
  const std::vector<RankingInformation::Link>& siteLinks
) {
  return siteToShapeVertexMap(stereopermutation, canonicalSites, siteLinks).invert();
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
