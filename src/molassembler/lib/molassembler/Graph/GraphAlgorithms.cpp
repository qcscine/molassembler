/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Graph/GraphAlgorithms.h"

#include "boost/graph/visitors.hpp"
#include "boost/graph/named_function_params.hpp"
#include "boost/graph/breadth_first_search.hpp"
#include "boost/graph/connected_components.hpp"
#include "boost/graph/biconnected_components.hpp"
#include "boost/range/iterator_range_core.hpp"
#include "boost/range/combine.hpp"
#include "temple/Adaptors/AllPairs.h"
#include "temple/Functional.h"
#include "temple/UnorderedSetAlgorithms.h"
#include "temple/TinySet.h"

#include "molassembler/AtomStereopermutator.h"
#include "molassembler/Cycles.h"
#include "molassembler/Modeling/AtomInfo.h"
#include "molassembler/RankingInformation.h"

namespace Scine {

namespace molassembler {

namespace graph_algorithms {

std::vector<LinkInformation> siteLinks(
  const PrivateGraph& graph,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB
) {
  // Make sure the stereopermutators called with are bonded in the first place
  assert(
    graph.edgeOption(
      stereopermutatorA.centralIndex(),
      stereopermutatorB.centralIndex()
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
      auto findEdgeIter = temple::find_if(
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

  std::vector<LinkInformation> links;
  std::map<
    std::pair<unsigned, unsigned>,
    unsigned
  > siteIndicesToLinksPositionMap;

  for(
    auto cycleOuterEdges : boost::make_iterator_range(
      graph.cycles().containing(
        BondIndex {
          stereopermutatorA.centralIndex(),
          stereopermutatorB.centralIndex()
        }
      )
    )
  ) {
    // Figure out which substituent of A and B is part of the cycle
    AtomIndex aAdjacent, bAdjacent;
    std::tie(aAdjacent, bAdjacent) = findNeighboringAtoms(
      stereopermutatorA.centralIndex(),
      stereopermutatorB.centralIndex(),
      cycleOuterEdges
    );

    std::pair<unsigned, unsigned> siteIndices {
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
      LinkInformation newLink;
      newLink.indexPair = siteIndices;
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

  temple::inplace::sort(links);
  return links;
}

std::vector<LinkInformation> siteLinks(
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
  std::unordered_map<AtomIndex, unsigned> indexToSiteMap;
  for(unsigned i = 0; i < sites.size(); ++i) {
    for(const AtomIndex siteAtomIndex : sites.at(i)) {
      indexToSiteMap.emplace(siteAtomIndex, i);
    }
  }

  std::vector<AtomIndex> sourceAdjacents;
  sourceAdjacents.reserve(graph.degree(source));
  for(const auto& siteAtomIndices : sites) {
    for(const AtomIndex adjacent : siteAtomIndices) {
      if(!temple::makeContainsPredicate(excludeAdjacents)(adjacent)) {
        sourceAdjacents.push_back(adjacent);
      }
    }
  }

  std::vector<LinkInformation> links;
  std::map<
    std::pair<unsigned, unsigned>,
    unsigned // Index of LinkInformation in links
  > siteIndicesToLinksPositionMap;

  temple::forEach(
    temple::adaptors::allPairs(sourceAdjacents),
    [&](const AtomIndex a, const AtomIndex b) {
      for(
        auto cycleOuterEdges :
        boost::make_iterator_range(
          graph.etaPreservedCycles().containing(
            std::vector<BondIndex> {
              BondIndex {source, a},
              BondIndex {source, b}
            }
          )
        )
      ) {
        assert(
          [&]() {
            auto containsPredicate = temple::makeContainsPredicate(cycleOuterEdges);
            return containsPredicate(BondIndex {source, a}) && containsPredicate(BondIndex {source, b});
          }()
        );

        const unsigned siteOfA = indexToSiteMap.at(a);
        const unsigned siteOfB = indexToSiteMap.at(b);

        if(siteOfA == siteOfB) {
          // If the cycle adjacents are from the same ligand, ignore this cycle
          continue;
        }

        const std::pair<unsigned, unsigned> siteIndices = std::minmax(siteOfA, siteOfB);

        const auto mapFindIter = siteIndicesToLinksPositionMap.find(siteIndices);
        if(
          mapFindIter == std::end(siteIndicesToLinksPositionMap)
          || cycleOuterEdges.size() < links.at(mapFindIter->second).cycleSequence.size()
        ) {
          auto newLink = LinkInformation {
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
  temple::inplace::sort(links);

  return links;
}

namespace detail {

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
   */
  return (
    siteAtoms.size() > 1
    && !( // Exclude the same-type triangle
      siteAtoms.size() == 2
      // The number of non-main-group elements is more than 1
      && temple::accumulate(
        siteAtoms,
        0u,
        [&](const unsigned carry, const AtomIndex adjacent) -> unsigned {
          if(!atom_info::isMainGroupElement(graph.elementType(adjacent))) {
            return carry + 1;
          }

          return carry;
        }
      ) > 1
    )
  );
}

void findSites(
  const PrivateGraph& graph,
  const AtomIndex centralIndex,
  const std::function<void(const std::vector<AtomIndex>&)>& callback
) {
  const unsigned A = graph.degree(centralIndex);
  temple::TinySet<PrivateGraph::Vertex> centralAdjacents;
  centralAdjacents.reserve(A);

  for(
    const PrivateGraph::Vertex adjacent :
    boost::make_iterator_range(graph.adjacents(centralIndex))
  ) {
    centralAdjacents.insert(adjacent);
  }

  std::vector<bool> skipList (A, false);

  temple::TinySet<AtomIndex> site;

  std::function<void(const PrivateGraph::Vertex)> recursiveDiscover
  = [&](const PrivateGraph::Vertex seed) {
    for(
      const PrivateGraph::Vertex adjacent :
      boost::make_iterator_range(graph.adjacents(seed))
    ) {
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

} // namespace detail

std::vector<
  std::vector<AtomIndex>
> ligandSiteGroups(
  const PrivateGraph& graph,
  AtomIndex centralIndex,
  const std::vector<AtomIndex>& excludeAdjacents
) {
  if(atom_info::isMainGroupElement(graph.elementType(centralIndex))) {
    std::vector<
      std::vector<AtomIndex>
    > adjacents;

    adjacents.reserve(graph.degree(centralIndex));

    for(
      const AtomIndex centralAdjacent :
      boost::make_iterator_range(graph.adjacents(centralIndex))
    ) {
      auto edge = graph.edge(centralAdjacent, centralIndex);

      if(graph.bondType(edge) == BondType::Eta) {
        /* A non-metal central index can have eta bonds, but they are not
         * considered in that atom's modeling, i.e. they do not form part of
         * their coordinate sphere.
         *
         * We determine whether the edge is mislabeled: Is the other vertex a
         * non-main-group element?
         */
        if(atom_info::isMainGroupElement(graph.elementType(centralAdjacent))) {
          std::string error = "Two main group elements are connected by an eta bond! ";
          error += std::to_string(centralAdjacent);
          error += " and ";
          error += std::to_string(centralIndex);

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
            centralIndex
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

  detail::findSites(
    graph,
    centralIndex,
    [&](const std::vector<AtomIndex>& ligand) -> void {
      // Make sure all bonds are marked properly
      if(detail::isHapticSite(ligand, graph)) {
        for(const auto& hapticIndex : ligand) {
          auto edge = graph.edge(centralIndex, hapticIndex);
          if(graph.bondType(edge) != BondType::Eta) {
            throw std::logic_error(
              "Haptic ligand constituting atom bound to non-main-group element via non-eta bond"
            );
          }
        }
      } else {
        for(const auto& nonHapticIndex : ligand) {
          auto edge = graph.edge(centralIndex, nonHapticIndex);
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
    temple::inplace::sort(ligand);
  }

  return groupedLigands;
}

void updateEtaBonds(PrivateGraph& graph) {
  const AtomIndex N = graph.N();
  for(AtomIndex centralIndex = 0; centralIndex < N; ++centralIndex) {
    // Skip any main group element types, none of these should be eta bonded
    if(atom_info::isMainGroupElement(graph.elementType(centralIndex))) {
      continue;
    }

    detail::findSites(
      graph,
      centralIndex,
      [&](const std::vector<AtomIndex>& ligand) -> void {
        if(detail::isHapticSite(ligand, graph)) {
          // Mark all bonds to the central atom as haptic bonds
          for(const auto& hapticIndex : ligand) {
            auto edge = graph.edge(centralIndex, hapticIndex);
            graph.bondType(edge) = BondType::Eta;
          }
        } else {
          // Mark all eta bonds to the central atom as single bonds
          for(const auto& hapticIndex : ligand) {
            auto edge = graph.edge(centralIndex, hapticIndex);
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
  assert(a < graph.N());

  std::vector<unsigned> distances (graph.N(), 0);

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

} // namespace graph_algorithms

} // namespace molassembler

} // namespace Scine
