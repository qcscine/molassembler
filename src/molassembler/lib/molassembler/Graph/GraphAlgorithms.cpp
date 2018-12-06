// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Graph/GraphAlgorithms.h"

#include "boost/graph/visitors.hpp"
#include "boost/graph/named_function_params.hpp"
#include "boost/graph/breadth_first_search.hpp"
#include "boost/graph/connected_components.hpp"
#include "boost/graph/biconnected_components.hpp"
#include "boost/range/iterator_range_core.hpp"
#include "boost/range/combine.hpp"
#include "temple/Functional.h"
#include "temple/UnorderedSetAlgorithms.h"
#include "temple/TinySet.h"

#include "molassembler/Modeling/AtomInfo.h"
#include "molassembler/Cycles.h"
#include "molassembler/RankingInformation.h"

namespace molassembler {

namespace GraphAlgorithms {

unsigned numConnectedComponents(const InnerGraph& graph) {
  std::vector<AtomIndex> component(graph.N());

  return boost::connected_components(graph.bgl(), &component[0]);
}

std::vector<LinkInformation> substituentLinks(
  const InnerGraph& graph,
  const Cycles& cycleData,
  const AtomIndex source,
  const std::vector<
    std::vector<AtomIndex>
  >& ligands,
  const std::vector<AtomIndex>& excludeAdjacents
) {
  // Early return: If the atom is terminal, there can be no links
  if(ligands.size() == 1) {
    return {};
  }

  /* General idea:
   *
   * In order to avoid a full O(N) BFS on every AtomStereopermutator candidate to
   * determine links, use cached Cycles instead to determine links.
   *
   * Iterate through every cycle and look for cycles containing exactly two of
   * the edge_descriptors corresponding to the edges from the source to an
   * active substituent.
   *
   * The two atoms adjacent to the central atom represented by those edges must
   * be from different ligands.
   */
  std::vector<LinkInformation> links;

  std::map<AtomIndex, unsigned> indexToLigandMap;
  for(unsigned i = 0; i < ligands.size(); ++i) {
    for(const auto& ligandIndex : ligands.at(i)) {
      indexToLigandMap.emplace(ligandIndex, i);
    }
  }

  temple::TinySet<InnerGraph::Edge> sourceAdjacentEdges;
  for(const auto& ligand : ligands) {
    for(const AtomIndex ligandConstitutingIndex : ligand) {
      if(
        std::find(
          std::begin(excludeAdjacents),
          std::end(excludeAdjacents),
          ligandConstitutingIndex
        ) == std::end(excludeAdjacents)
      ) {
        sourceAdjacentEdges.insert(
          graph.edge(source, ligandConstitutingIndex)
        );
      }
    }
  }

  std::map<
    std::pair<unsigned, unsigned>,
    unsigned // Index of LinkInformation in links
  > linksMap;

  for(const auto cyclePtr : cycleData) {
    auto cycleEdges = temple::map(
      Cycles::edges(cyclePtr),
      [&graph](const BondIndex& bond) -> InnerGraph::Edge {
        return graph.edge(bond.first, bond.second);
      }
    );

    // For set_intersection to work, cycle edges need to be sorted too
    std::sort(
      std::begin(cycleEdges),
      std::end(cycleEdges)
    );

    std::vector<InnerGraph::Edge> intersection;

    // Set intersection copies intersected elements from the first range
    std::set_intersection(
      std::begin(sourceAdjacentEdges),
      std::end(sourceAdjacentEdges),
      std::begin(cycleEdges),
      std::end(cycleEdges),
      std::back_inserter(intersection)
    );

    // Must match exactly two edges
    if(intersection.size() == 2) {
      /* Intersected edges come from sourceAdjacentEdges, which are from an
       * out_edges call, where it is guaranteed that source is the vertex
       * out_edges was called with. The atom indices of interest are therefore
       * the targets of these edges.
       */
      const AtomIndex a = graph.target(intersection.front());
      const AtomIndex b = graph.target(intersection.back());

      const unsigned ligandOfA = indexToLigandMap.at(a);
      const unsigned ligandOfB = indexToLigandMap.at(b);

      if(ligandOfA == ligandOfB) {
        // If the cycle adjacents are from the same ligand, ignore this cycle
        continue;
      }

      auto indexPair = std::pair<unsigned, unsigned> {
        std::min(ligandOfA, ligandOfB),
        std::max(ligandOfA, ligandOfB),
      };

      // Find which adjacents are overlapping
      if(
        linksMap.count(indexPair) == 0
        || (
          linksMap.count(indexPair) > 0
          && cycleEdges.size() < links.at(linksMap.at(indexPair)).cycleSequence.size()
        )
      ) {
        auto newLink = LinkInformation {
          indexPair,
          makeRingIndexSequence(Cycles::edgeVertices(cyclePtr)),
          source
        };

        /* Track the current best link for a given ligand index pair and
         * improve it if a better one turns up
         */
        if(linksMap.count(indexPair) == 0) {
          links.push_back(std::move(newLink));
          linksMap.emplace(
            indexPair,
            links.size() - 1
          );
        } else {
          links.at(linksMap.at(indexPair)) = std::move(newLink);
        }
      }
    }
  }

  // Sort the links before passing them out in order to ease comparisons
  temple::inplace::sort(links);

  return links;
}

namespace detail {

bool isHapticLigand(
  const std::vector<AtomIndex>& ligand,
  const InnerGraph& graph
) {
  /* A ligand is a haptic ligand if:
   * - The number of co-bonded atoms constituting direct bonds is more than one
   * - It is not a same-type triangle, i.e. we have detected
   *
   *     M — B
   *      \ /   where M, A and B are all non-main-group elements
   *       A
   */
  return (
    ligand.size() > 1
    && !( // Exclude the same-type triangle
      ligand.size() == 2
      // The number of non-main-group elements is more than 1
      && temple::accumulate(
        ligand,
        0u,
        [&](const unsigned carry, const AtomIndex adjacent) -> unsigned {
          if(!AtomInfo::isMainGroupElement(graph.elementType(adjacent))) {
            return carry + 1;
          }

          return carry;
        }
      ) > 1
    )
  );
}

void findLigands(
  const InnerGraph& graph,
  const AtomIndex centralIndex,
  const std::function<void(const std::vector<AtomIndex>&)>& callback
) {
  unsigned A = graph.degree(centralIndex);
  temple::TinySet<InnerGraph::Vertex> centralAdjacents;
  centralAdjacents.reserve(A);

  for(
    const InnerGraph::Vertex adjacent :
    boost::make_iterator_range(graph.adjacents(centralIndex))
  ) {
    centralAdjacents.insert(adjacent);
  }

  std::vector<bool> skipList (A, false);

  temple::TinySet<AtomIndex> ligand;

  std::function<void(const InnerGraph::Vertex)> recursiveDiscover
  = [&](const InnerGraph::Vertex seed) {
    for(
      const InnerGraph::Vertex adjacent :
      boost::make_iterator_range(graph.adjacents(seed))
    ) {
      if(centralAdjacents.count(adjacent) > 0 && ligand.count(adjacent) == 0) {
        // *iter is shared adjacent of center and seed and not yet discovered
        ligand.insert(adjacent);
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
    ligand = {adjacent};
    recursiveDiscover(adjacent);

    callback(ligand.data);
  }
}

} // namespace detail

std::vector<
  std::vector<AtomIndex>
> ligandSiteGroups(
  const InnerGraph& graph,
  AtomIndex centralIndex,
  const std::vector<AtomIndex>& excludeAdjacents
) {
  // A non-metal central index cannot have any eta bonds
  if(AtomInfo::isMainGroupElement(graph.elementType(centralIndex))) {
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
        // Determine if the edge is mislabeled
        // is the other vertex a non-main-group element?

        if(AtomInfo::isMainGroupElement(graph.elementType(centralAdjacent))) {
          std::string error = "Two main group elements are connected by an eta bond! ";
          error += std::to_string(centralAdjacent);
          error += " and ";
          error += std::to_string(centralIndex);

          throw std::logic_error(error);
        }
        /* If the central index is a main-group element and connected to a
         * non-main-group element via an eta bond, this adjacency is not
         * recorded for the determination of ligand site groups. The symmetry
         * modelling at that center should not include individual eta bond
         * contributions.
         */
      } else {
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

  detail::findLigands(
    graph,
    centralIndex,
    [&](const std::vector<AtomIndex>& ligand) -> void {
      // Make sure all bonds are marked properly
      if(detail::isHapticLigand(ligand, graph)) {
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

  return groupedLigands;
}

void findAndSetEtaBonds(InnerGraph& graph) {
  const AtomIndex N = graph.N();
  for(AtomIndex centralIndex = 0; centralIndex < N; ++centralIndex) {
    // Skip any main group element types, none of these should be eta bonded
    if(AtomInfo::isMainGroupElement(graph.elementType(centralIndex))) {
      continue;
    }

    detail::findLigands(
      graph,
      centralIndex,
      [&](const std::vector<AtomIndex>& ligand) -> void {
        if(detail::isHapticLigand(ligand, graph)) {
          // Mark all bonds to the central atom as haptic bonds
          for(const auto& hapticIndex : ligand) {
            auto edge = graph.edge(centralIndex, hapticIndex);
            graph.bondType(edge) = BondType::Eta;
          }
        }
      }
    );
  }
}

std::vector<unsigned> distance(AtomIndex a, const InnerGraph& graph) {
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

} // namespace GraphAlgorithms

} // namespace molassembler
