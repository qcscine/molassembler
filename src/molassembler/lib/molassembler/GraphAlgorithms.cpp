#include "boost/graph/connected_components.hpp"
#include "boost/graph/biconnected_components.hpp"
#include "boost/range/combine.hpp"

#include "Delib/ElementInfo.h"

#include "temple/Containers.h"
#include "temple/UnorderedSets.h"
#include "temple/TinySet.h"

#include "AtomInfo.h"
#include "Cycles.h"
#include "GraphAlgorithms.h"

#include <iostream>

namespace molassembler {

namespace GraphAlgorithms {

LinkInformation::LinkInformation() = default;

LinkInformation::LinkInformation(
  std::pair<unsigned, unsigned> ligandIndices,
  std::vector<AtomIndexType> sequence,
  const AtomIndexType source
) {
  /* Fix degrees of freedom of the underlying information so we can
   * efficiently implement operator <. indexPair can be an ordered pair:
   */
  indexPair = std::move(ligandIndices);
  if(indexPair.first > indexPair.second) {
    std::swap(indexPair.first, indexPair.second);
  }

  // The cycle sequence should be centralized on the source vertex
  cycleSequence = centralizeRingIndexSequence(std::move(sequence), source);

  /* After centralization, the source vertex is first and last. We need to
   * fix the remaining degree of freedom, which is if that the cycle sequence
   * in between can be reversed. We choose to fix it by making it ascending
   * if there are at least four vertices in the sequence between the second
   * and second-to-last vertices
   */

  if(
    cycleSequence.size() >= 4
    && cycleSequence.at(1) > cycleSequence.at(cycleSequence.size() - 2)
  ) {
    // Reverse is [first, last), and cycleSequence is a vector, so:
    std::reverse(
      std::begin(cycleSequence) + 1,
      std::end(cycleSequence) - 1
    );
  }
}

bool LinkInformation::operator == (const LinkInformation& other) const {
  return (
    indexPair == other.indexPair
    && cycleSequence == other.cycleSequence
  );
}

bool LinkInformation::operator != (const LinkInformation& other) const {
  return !(*this == other);
}

bool LinkInformation::operator < (const LinkInformation& other) const {
  return std::tie(indexPair, cycleSequence) < std::tie(other.indexPair, other.cycleSequence);
}

unsigned numConnectedComponents(const GraphType& graph) {
  std::vector<AtomIndexType> component(boost::num_vertices(graph));

  return boost::connected_components(graph, &component[0]);
}

std::vector<LinkInformation> substituentLinks(
  const GraphType& graph,
  const Cycles& cycleData,
  const AtomIndexType source,
  const std::vector<
    std::vector<AtomIndexType>
  >& ligands,
  const std::set<AtomIndexType>& excludeAdjacents
) {
  /* General idea:
   *
   * In order to avoid a full O(N) BFS on every AtomStereocenter candidate to
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

  std::map<AtomIndexType, unsigned> indexToLigandMap;
  for(unsigned i = 0; i < ligands.size(); ++i) {
    for(const auto& ligandIndex : ligands.at(i)) {
      indexToLigandMap.emplace(
        ligandIndex,
        i
      );
    }
  }

  temple::TinySet<GraphType::edge_descriptor> sourceAdjacentEdges;
  for(const auto& ligand : ligands) {
    for(const AtomIndexType ligandConstitutingIndex : ligand) {
      if(excludeAdjacents.count(ligandConstitutingIndex) == 0) {
        sourceAdjacentEdges.insert(
          boost::edge(source, ligandConstitutingIndex, graph).first
        );
      }
    }
  }

  std::map<
    std::pair<unsigned, unsigned>,
    unsigned // Index of LinkInformation in links
  > linksMap;

  auto getAdjacentOfEdge = [&](const auto& edgeDescriptor) -> AtomIndexType {
    auto edgeSource = boost::source(edgeDescriptor, graph);

    if(source == edgeSource) {
      return boost::target(edgeDescriptor, graph);
    }

    return edgeSource;
  };

  for(const auto cyclePtr : cycleData) {
    auto cycleEdges = Cycles::edges(cyclePtr, graph);

    std::vector<GraphType::edge_descriptor> intersection;

    std::set_intersection(
      sourceAdjacentEdges.begin(),
      sourceAdjacentEdges.end(),
      cycleEdges.begin(),
      cycleEdges.end(),
      std::back_inserter(intersection)
    );

    // Must match exactly two edges
    if(intersection.size() == 2) {
      AtomIndexType a = getAdjacentOfEdge(intersection.front());
      AtomIndexType b = getAdjacentOfEdge(intersection.back());

      unsigned ligandOfA = indexToLigandMap.at(a);
      unsigned ligandOfB = indexToLigandMap.at(b);

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
  temple::sort(links);

  return links;
}

namespace detail {

bool isHapticLigand(
  const std::vector<AtomIndexType>& ligand,
  const GraphType& graph
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
        [&](const unsigned carry, const AtomIndexType adjacent) -> unsigned {
          if(!AtomInfo::isMainGroupElement(graph[adjacent].elementType)) {
            return carry + 1;
          }

          return carry;
        }
      ) > 1
    )
  );
}

void findLigands(
  const GraphType& graph,
  const AtomIndexType centralIndex,
  const std::function<void(const std::vector<AtomIndexType>&)>& callback
) {
  unsigned A = boost::out_degree(centralIndex, graph);
  temple::TinySet<AtomIndexType> centralAdjacents;

  GraphType::adjacency_iterator iter, end;
  std::tie(iter, end) = boost::adjacent_vertices(centralIndex, graph);

  for(; iter != end; ++iter) {
    centralAdjacents.insert(*iter);
  }

  std::vector<bool> skipList (A, false);

  temple::TinySet<AtomIndexType> ligand;

  std::function<void(const AtomIndexType&)> recursiveDiscover = [&](const AtomIndexType& seed) {
    GraphType::adjacency_iterator iter, end;
    std::tie(iter, end) = boost::adjacent_vertices(seed, graph);
    for(; iter != end; ++iter) {
      if(centralAdjacents.count(*iter) > 0 && ligand.count(*iter) == 0) {
        // *iter is shared adjacent of center and seed and not yet discovered
        ligand.insert(*iter);
        recursiveDiscover(*iter);
      }
    }

    // Make sure the algorithm skips this seed index to avoid duplicate work
    skipList.at(
      centralAdjacents.find(seed) - centralAdjacents.begin()
    ) = true;
  };

  for(unsigned i = 0; i < A; ++i) {
    if(skipList.at(i)) {
      continue;
    }

    const AtomIndexType& adjacent = centralAdjacents.at(i);
    ligand = {adjacent};
    recursiveDiscover(adjacent);

    callback(ligand.data);
  }
}

} // namespace detail

std::vector<
  std::vector<AtomIndexType>
> ligandSiteGroups(
  const GraphType& graph,
  AtomIndexType centralIndex,
  const std::set<AtomIndexType>& excludeAdjacents
) {
  // A non-metal central index cannot have any eta bonds
  if(AtomInfo::isMainGroupElement(graph[centralIndex].elementType)) {
    std::vector<
      std::vector<AtomIndexType>
    > adjacents;

    adjacents.reserve(boost::out_degree(centralIndex, graph));

    for(
      const AtomIndexType centralAdjacent :
      RangeForTemporary<GraphType::adjacency_iterator>(
        boost::adjacent_vertices(centralIndex, graph)
      )
    ) {
      auto edgeFoundPair = boost::edge(centralAdjacent, centralIndex, graph);

      if(graph[edgeFoundPair.first].bondType == BondType::Eta) {
        // Determine if the edge is mislabeled
        // is the other vertex a non-main-group element?

        if(AtomInfo::isMainGroupElement(graph[centralAdjacent].elementType)) {
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
        if(excludeAdjacents.count(centralAdjacent) == 0) {
          adjacents.push_back(
            std::vector<AtomIndexType> {centralAdjacent}
          );
        }
      }
    }

    return adjacents;
  }

  std::vector<
    std::vector<AtomIndexType>
  > groupedLigands;

  detail::findLigands(
    graph,
    centralIndex,
    [&](const std::vector<AtomIndexType>& ligand) -> void {
      // Make sure all bonds are marked properly
      if(detail::isHapticLigand(ligand, graph)) {
        for(const auto& hapticIndex : ligand) {
          auto edgePair = boost::edge(centralIndex, hapticIndex, graph);
          if(graph[edgePair.first].bondType != BondType::Eta) {
            throw std::logic_error(
              "Haptic ligand constituting atom bound to non-main-group element via non-eta bond"
            );
          }
        }
      } else {
        for(const auto& nonHapticIndex : ligand) {
          auto edgePair = boost::edge(centralIndex, nonHapticIndex, graph);
          if(graph[edgePair.first].bondType == BondType::Eta) {
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
          && excludeAdjacents.count(ligand.front()) == 1
        )
      ) {
        // Add the ligand to the set
        groupedLigands.emplace_back(ligand.begin(), ligand.end());
      }
    }
  );

  return groupedLigands;
}

void findAndSetEtaBonds(GraphType& graph) {
  AtomIndexType N = boost::num_vertices(graph);
  for(AtomIndexType centralIndex = 0; centralIndex < N; ++centralIndex) {
    // Skip any main group element types, none of these should be eta bonded
    if(AtomInfo::isMainGroupElement(graph[centralIndex].elementType)) {
      break;
    }

    detail::findLigands(
      graph,
      centralIndex,
      [&](const std::vector<AtomIndexType>& ligand) -> void {
        if(detail::isHapticLigand(ligand, graph)) {
          // Mark all bonds to the central atom as haptic bonds
          for(const auto& hapticIndex : ligand) {
            auto edgePair = boost::edge(centralIndex, hapticIndex, graph);
            graph[edgePair.first].bondType = BondType::Eta;
          }
        }
      }
    );
  }
}

RemovalSafetyData getRemovalSafetyData(const GraphType& graph) {
  RemovalSafetyData safetyData;

  std::vector<AtomIndexType> articulationVertices;

  using ComponentMapBase = std::map<
    GraphType::edge_descriptor,
    GraphType::edges_size_type
  >;

  ComponentMapBase componentMapData;
  boost::associative_property_map<ComponentMapBase> componentMap(componentMapData);
  unsigned numComponents;

  // Calculate the biconnected components and articulation vertices
  std::tie(numComponents, std::ignore) = boost::biconnected_components(
    graph,
    componentMap,
    std::back_inserter(articulationVertices)
  );

  // Copy articulation vertices to the set
  for(const auto& vertex : articulationVertices) {
    safetyData.articulationVertices.insert(vertex);
  }

  /* Work out from the biconnected components which edges are bridges: If the
   * biconnected component encompasses only a single edge, it is a bridge
   */
  std::vector<
    std::set<EdgeIndexType>
  > componentSets (numComponents);

  for(const auto& mapIterPair : componentMapData) {
    const auto& edgeIndex = mapIterPair.first;
    const auto& componentIndex = mapIterPair.second;

    componentSets[componentIndex].insert(edgeIndex);
  }

  for(const auto& componentSet : componentSets) {
    if(componentSet.size() == 1) {
      safetyData.bridges.insert(
        *componentSet.begin()
      );
    }
  }

  return safetyData;
}

} // namespace GraphAlgorithms

} // namespace molassembler
