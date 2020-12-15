/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Graph/GraphAlgorithms.h"

#include "boost/graph/visitors.hpp"
#include "boost/graph/named_function_params.hpp"
#include "boost/graph/breadth_first_search.hpp"
#include "boost/graph/connected_components.hpp"
#include "boost/graph/biconnected_components.hpp"
#include "boost/range/combine.hpp"
#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/UnorderedSetAlgorithms.h"
#include "Molassembler/Temple/TinySet.h"

#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/Cycles.h"
#include "Molassembler/Modeling/AtomInfo.h"
#include "Molassembler/RankingInformation.h"

namespace Scine {
namespace Molassembler {
namespace GraphAlgorithms {

std::vector<RankingInformation::Link> siteLinks(
  const PrivateGraph& graph,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB
) {
  // Make sure the stereopermutators called with are bonded in the first place
  assert(
    graph.edgeOption(
      stereopermutatorA.placement(),
      stereopermutatorB.placement()
    )
  );

  // If either permutator is terminal, there can be no links
  if(
    stereopermutatorA.getRanking().sites.size() == 1
    || stereopermutatorB.getRanking().sites.size() == 1
  ) {
    return {};
  }

  auto findNeighboringAtoms = [](
    const AtomIndex a,
    const AtomIndex b,
    const std::vector<BondIndex>& cycleBonds
  ) -> std::pair<AtomIndex, AtomIndex> {
    auto findNeighbor = [](
      const AtomIndex focalAtom,
      const AtomIndex avoidAtom,
      const std::vector<BondIndex>& cycle
    ) -> AtomIndex {
      auto findEdgeIter = Temple::find_if(
        cycle,
        [&](const BondIndex& edge) -> bool {
          return (edge.contains(focalAtom) && !edge.contains(avoidAtom));
        }
      );

      assert(findEdgeIter != std::end(cycle));

      return (
        findEdgeIter->first == focalAtom
        ? findEdgeIter->second
        : findEdgeIter->first
      );
    };

    return {
      findNeighbor(a, b, cycleBonds),
      findNeighbor(b, a, cycleBonds)
    };
  };

  std::vector<RankingInformation::Link> links;
  std::map<
    std::pair<SiteIndex, SiteIndex>,
    unsigned
  > siteIndicesToLinksPositionMap;

  for(
    auto cycleOuterEdges :
    graph.cycles().containing(
      BondIndex {
        stereopermutatorA.placement(),
        stereopermutatorB.placement()
      }
    )
  ) {
    // Figure out which substituent of A and B is part of the cycle
    AtomIndex aAdjacent {};
    AtomIndex bAdjacent {};
    std::tie(aAdjacent, bAdjacent) = findNeighboringAtoms(
      stereopermutatorA.placement(),
      stereopermutatorB.placement(),
      cycleOuterEdges
    );

    std::pair<SiteIndex, SiteIndex> siteIndices {
      stereopermutatorA.getRanking().getSiteIndexOf(aAdjacent),
      stereopermutatorB.getRanking().getSiteIndexOf(bAdjacent)
    };

    const auto mapFindIter = siteIndicesToLinksPositionMap.find(siteIndices);
    if(
      // No entry for site indices
      mapFindIter == std::end(siteIndicesToLinksPositionMap)
      // Or new cycle is smaller than one already found
      || cycleOuterEdges.size() < links.at(mapFindIter->second).cycleSequence.size()
    ) {
      // No invariant establishing for these LinkInformations
      RankingInformation::Link newLink;
      newLink.sites = siteIndices;
      newLink.cycleSequence = makeRingIndexSequence(std::move(cycleOuterEdges));

      // Add or update existing
      if(mapFindIter == std::end(siteIndicesToLinksPositionMap)) {
        links.push_back(std::move(newLink));
        siteIndicesToLinksPositionMap.emplace(
          siteIndices,
          links.size() - 1
        );
      } else {
        links.at(mapFindIter->second) = std::move(newLink);
      }
    }
  }

  Temple::sort(links);
  return links;
}

std::vector<RankingInformation::Link> siteLinks(
  const PrivateGraph& graph,
  const AtomIndex source,
  const std::vector<
    std::vector<AtomIndex>
  >& sites,
  const std::vector<AtomIndex>& excludeAdjacents
) {
  // Early return: If the atom is terminal, there can be no links
  if(sites.size() == 1) {
    return {};
  }

  /* General idea:
   *
   * In order to avoid a full O(N) BFS on every AtomStereopermutator candidate
   * to determine links, use Cycles instead.
   *
   * Iterate through every cycle and look for cycles containing two of source's
   * immediate adjacents.
   *
   * The two atoms adjacent to the central atom represented by those edges must
   * be from different sites.
   */
  std::unordered_map<AtomIndex, SiteIndex> indexToSiteMap;
  for(unsigned i = 0; i < sites.size(); ++i) {
    for(const AtomIndex siteAtomIndex : sites.at(i)) {
      indexToSiteMap.emplace(siteAtomIndex, i);
    }
  }

  std::vector<AtomIndex> sourceAdjacents;
  sourceAdjacents.reserve(graph.degree(source));
  for(const auto& siteAtomIndices : sites) {
    for(const AtomIndex adjacent : siteAtomIndices) {
      if(!Temple::makeContainsPredicate(excludeAdjacents)(adjacent)) {
        sourceAdjacents.push_back(adjacent);
      }
    }
  }

  std::vector<RankingInformation::Link> links;
  std::map<
    std::pair<SiteIndex, SiteIndex>,
    unsigned // Index of LinkInformation in links
  > siteIndicesToLinksPositionMap;

  Temple::forEach(
    Temple::Adaptors::allPairs(sourceAdjacents),
    [&](const AtomIndex a, const AtomIndex b) {
      for(
        auto cycleOuterEdges :
        graph.etaPreservedCycles().containing(
          std::vector<BondIndex> {
            BondIndex {source, a},
            BondIndex {source, b}
          }
        )
      ) {
        assert(
          [&]() {
            auto containsPredicate = Temple::makeContainsPredicate(cycleOuterEdges);
            return containsPredicate(BondIndex {source, a}) && containsPredicate(BondIndex {source, b});
          }()
        );

        const SiteIndex siteOfA = indexToSiteMap.at(a);
        const SiteIndex siteOfB = indexToSiteMap.at(b);

        if(siteOfA == siteOfB) {
          // If the cycle adjacents are from the same ligand, ignore this cycle
          continue;
        }

        const std::pair<SiteIndex, SiteIndex> siteIndices = std::minmax(siteOfA, siteOfB);

        const auto mapFindIter = siteIndicesToLinksPositionMap.find(siteIndices);
        if(
          mapFindIter == std::end(siteIndicesToLinksPositionMap)
          || cycleOuterEdges.size() < links.at(mapFindIter->second).cycleSequence.size()
        ) {
          auto newLink = RankingInformation::Link {
            siteIndices,
            makeRingIndexSequence(std::move(cycleOuterEdges)),
            source
          };

          /* Track the current best link for a given ligand index pair and
           * improve it if a better one turns up
           */
          if(mapFindIter == std::end(siteIndicesToLinksPositionMap)) {
            links.push_back(std::move(newLink));
            siteIndicesToLinksPositionMap.emplace(
              siteIndices,
              links.size() - 1
            );
          } else {
            links.at(mapFindIter->second) = std::move(newLink);
          }
        }
      }
    }
  );

  // Sort the links before passing them out in order to ease comparisons
  Temple::sort(links);
  return links;
}

namespace {

bool isHapticSite(
  const std::vector<AtomIndex>& siteAtoms,
  const PrivateGraph& graph
) {
  /* A site is haptic if:
   * - The number of co-bonded atoms constituting direct bonds is more than one
   * - And it is not a same-type triangle, i.e. we have detected
   *
   *     M — B
   *      \ /   where M, A and B are all non-main-group elements
   *       A
   *
   * The number of non-main-group elements is at most one
   */
  return (
    siteAtoms.size() > 1
    && Temple::accumulate(
      siteAtoms,
      0U,
      [&](const unsigned carry, const AtomIndex adjacent) -> unsigned {
        if(!AtomInfo::isMainGroupElement(graph.elementType(adjacent))) {
          return carry + 1;
        }

        return carry;
      }
    ) <= 1
  );
}

void findSites(
  const PrivateGraph& graph,
  const AtomIndex placement,
  const std::function<void(const std::vector<AtomIndex>&)>& callback
) {
  const unsigned A = graph.degree(placement);
  Temple::TinySet<PrivateGraph::Vertex> centralAdjacents;
  centralAdjacents.reserve(A);

  for(const PrivateGraph::Vertex adjacent : graph.adjacents(placement)) {
    centralAdjacents.insert(adjacent);
  }

  std::vector<bool> skipList (A, false);

  Temple::TinySet<AtomIndex> site;

  std::function<void(const PrivateGraph::Vertex)> recursiveDiscover
  = [&](const PrivateGraph::Vertex seed) {
    for(const PrivateGraph::Vertex adjacent : graph.adjacents(seed)) {
      if(centralAdjacents.count(adjacent) > 0 && site.count(adjacent) == 0) {
        // *iter is shared adjacent of center and seed and not yet discovered
        site.insert(adjacent);
        recursiveDiscover(adjacent);
      }
    }

    // Make sure the algorithm skips this seed index to avoid duplicate work
    skipList.at(
      centralAdjacents.find(seed) - centralAdjacents.begin()
    ) = true;
  };

  // Go through all adjacent indices
  for(unsigned i = 0; i < A; ++i) {
    if(skipList.at(i)) {
      continue;
    }

    const AtomIndex& adjacent = centralAdjacents.at(i);
    site = {adjacent};
    recursiveDiscover(adjacent);

    callback(site.data);
  }
}

} // namespace

std::vector<
  std::vector<AtomIndex>
> sites(
  const PrivateGraph& graph,
  AtomIndex placement,
  const std::vector<AtomIndex>& excludeAdjacents
) {
  if(AtomInfo::isMainGroupElement(graph.elementType(placement))) {
    std::vector<
      std::vector<AtomIndex>
    > adjacents;

    adjacents.reserve(graph.degree(placement));

    for(const AtomIndex centralAdjacent : graph.adjacents(placement)) {
      auto edge = graph.edge(centralAdjacent, placement);

      if(graph.bondType(edge) == BondType::Eta) {
        /* A non-metal central index can have eta bonds, but they are not
         * considered in that atom's modeling, i.e. they do not form part of
         * their coordinate sphere.
         *
         * We determine whether the edge is mislabeled: Is the other vertex a
         * non-main-group element?
         */
        if(AtomInfo::isMainGroupElement(graph.elementType(centralAdjacent))) {
          std::string error = "Two main group elements are connected by an eta bond! ";
          error += std::to_string(centralAdjacent);
          error += " and ";
          error += std::to_string(placement);

          throw std::logic_error(error);
        }
        /* If the central index is a main-group element and connected to a
         * non-main-group element via an eta bond, this adjacency is not
         * recorded for the determination of ligand site groups. The shape
         * modelling at that center should not include individual eta bond
         * contributions.
         */
      } else {
        // This is a regular edge
        if(
          std::find(
            std::begin(excludeAdjacents),
            std::end(excludeAdjacents),
            placement
          ) == std::end(excludeAdjacents)
        ) {
          adjacents.push_back(
            std::vector<AtomIndex> {centralAdjacent}
          );
        }
      }
    }

    return adjacents;
  }

  std::vector<
    std::vector<AtomIndex>
  > groupedLigands;

  findSites(
    graph,
    placement,
    [&](const std::vector<AtomIndex>& ligand) -> void {
      // Make sure all bonds are marked properly
      if(isHapticSite(ligand, graph)) {
        for(const auto& hapticIndex : ligand) {
          auto edge = graph.edge(placement, hapticIndex);
          if(graph.bondType(edge) != BondType::Eta) {
            throw std::logic_error(
              "Haptic ligand constituting atom bound to non-main-group element via non-eta bond"
            );
          }
        }
      } else {
        for(const auto& nonHapticIndex : ligand) {
          auto edge = graph.edge(placement, nonHapticIndex);
          if(graph.bondType(edge) == BondType::Eta) {
            throw std::logic_error(
              "Non-haptic ligand bound to non-main-group element via eta bond"
            );
          }
        }
      }

      /* We need to exclude this ligand if it is size one and consists of an
       * excluded adjacent
       */
      if(
        !(
          ligand.size() == 1
          && std::find(
            std::begin(excludeAdjacents),
            std::end(excludeAdjacents),
            ligand.front()
          ) != std::end(excludeAdjacents)
        )
      ) {
        // Add the ligand to the set
        groupedLigands.emplace_back(ligand.begin(), ligand.end());
      }
    }
  );

  for(auto& ligand : groupedLigands) {
    Temple::sort(ligand);
  }

  return groupedLigands;
}

void updateEtaBonds(PrivateGraph& graph) {
  const AtomIndex N = graph.V();
  for(AtomIndex placement = 0; placement < N; ++placement) {
    // Skip any main group element types, none of these should be eta bonded
    if(AtomInfo::isMainGroupElement(graph.elementType(placement))) {
      continue;
    }

    findSites(
      graph,
      placement,
      [&](const std::vector<AtomIndex>& ligand) -> void {
        if(isHapticSite(ligand, graph)) {
          // Mark all bonds to the central atom as haptic bonds
          for(const auto& hapticIndex : ligand) {
            auto edge = graph.edge(placement, hapticIndex);
            graph.bondType(edge) = BondType::Eta;
          }
        } else {
          // Mark all eta bonds to the central atom as single bonds
          for(const auto& hapticIndex : ligand) {
            auto edge = graph.edge(placement, hapticIndex);
            BondType& bondType = graph.bondType(edge);
            if(bondType == BondType::Eta) {
              bondType = BondType::Single;
            }
          }
        }
      }
    );
  }
}

std::vector<unsigned> distance(AtomIndex a, const PrivateGraph& graph) {
  assert(a < graph.V());

  std::vector<unsigned> distances (graph.V(), 0);

  boost::breadth_first_search(
    graph.bgl(),
    a,
    boost::visitor(
      boost::make_bfs_visitor(
        boost::record_distances(&distances[0], boost::on_tree_edge {})
      )
    )
  );

  return distances;
}

std::vector<AtomIndex> shortestPaths(AtomIndex a, const PrivateGraph& graph) {
  assert(a < graph.V());

  std::vector<AtomIndex> predecessors(graph.V(), 0);

  boost::breadth_first_search(
    graph.bgl(),
    a,
    boost::visitor(
      boost::make_bfs_visitor(
        boost::record_predecessors(&predecessors[0], boost::on_tree_edge())
      )
    )
  );

  // Link the source vertex to itself
  predecessors.at(a) = a;
  return predecessors;
}

std::vector<AtomIndex> bfsUniqueDescendants(
  const AtomIndex source,
  const std::vector<AtomIndex>& descendants,
  const PrivateGraph& graph
) {
  assert(source < graph.V());
  constexpr AtomIndex unmarked = std::numeric_limits<AtomIndex>::max() - 1;
  constexpr AtomIndex split = std::numeric_limits<AtomIndex>::max();

  struct Visitor {
    using Graph = PrivateGraph::BglType;
    using Vertex = Graph::vertex_descriptor;
    using Edge = Graph::edge_descriptor;

    explicit Visitor(const AtomIndex a, std::vector<AtomIndex>& components, const Graph& g)
      : componentMembership(components)
    {
      depths.resize(boost::num_vertices(g), 0);
      depths.at(a) = 0;
    }

    void initialize_vertex(const Vertex /* v */, const Graph& /* g */) {}
    void discover_vertex(const Vertex /* v */, const Graph& /* g */) {}
    void examine_vertex(const Vertex /* v */, const Graph& /* g */) {}
    void examine_edge(const Edge /* e */, const Graph& /* g */) {}
    void tree_edge(const Edge e, const Graph& g) {
      const Vertex source = boost::source(e, g);
      const Vertex target = boost::target(e, g);
      auto& components = componentMembership.get();
      if(components.at(target) == unmarked) {
        components.at(target) = components.at(source);
      }
      depths.at(target) = depths.at(source) + 1;
    }
    void non_tree_edge(const PrivateGraph::Edge /* e */, const Graph& /* g */) {}
    void gray_target(const PrivateGraph::Edge e, const Graph& g) {
      auto& components = componentMembership.get();
      const Vertex source = boost::source(e, g);
      const Vertex target = boost::target(e, g);
      if(
        depths.at(target) > depths.at(source)
        && components.at(target) != components.at(source)
      ) {
        components.at(target) = split;
      }
    }
    void black_target(const PrivateGraph::Edge /* e */, const Graph& /* g */) {}
    void finish_vertex(const PrivateGraph::Vertex /* v */, const Graph& /* g */) {}

    std::vector<unsigned> depths;
    std::reference_wrapper<std::vector<AtomIndex>> componentMembership;
  };

  std::vector<AtomIndex> components(graph.V(), unmarked);
  components.at(source) = source;
  for(const AtomIndex descendant : descendants) {
    components.at(descendant) = descendant;
  }

  boost::breadth_first_search(
    graph.bgl(),
    source,
    boost::visitor(
      Visitor(source, components, graph.bgl())
    )
  );

  assert(Temple::all_of(components, [&](const AtomIndex grp) { return grp != unmarked; }));

  return components;
}

} // namespace GraphAlgorithms
} // namespace Molassembler
} // namespace Scine
