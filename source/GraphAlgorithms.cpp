#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>

#include "CycleData.h"
#include "GraphAlgorithms.h"

#include <iostream>

namespace MoleculeManip {

namespace GraphAlgorithms {

unsigned numConnectedComponents(const GraphType& graph) {
  std::vector<AtomIndexType> component(boost::num_vertices(graph));

  return boost::connected_components(graph, &component[0]);
}

std::vector<LinkInformation> substituentLinks(
  const GraphType& graph,
  const CycleData& cycleData,
  const AtomIndexType& source,
  const std::vector<AtomIndexType>& activeAdjacents
) {
  /* General idea:
   *
   * In order to avoid a full O(N) BFS on every CNStereocenter candidate to
   * determine links, use cached CycleData instead to determine links.
   *
   * Iterate through every cycle and look for cycles containing exactly two of
   * the edge_descriptors corresponding to the edges from the source to an
   * active substituent
   */
  std::vector<LinkInformation> links;

  temple::TinySet<GraphType::edge_descriptor> sourceAdjacentEdges;
  for(const auto& adjacentIndex : activeAdjacents) {
    sourceAdjacentEdges.insert(
      boost::edge(
        source,
        adjacentIndex,
        graph
      ).first
    );
  }

  std::map<
    std::pair<AtomIndexType, AtomIndexType>,
    unsigned
  > linksMap;

  auto getAdjacentOfEdge = [&](const auto& edgeDescriptor) -> AtomIndexType {
    auto edgeSource = boost::source(edgeDescriptor, graph);

    if(source == edgeSource) {
      return boost::target(edgeDescriptor, graph);
    }

    return edgeSource;
  };

  for(
    auto cycleIterator = cycleData.getCyclesIterator();
    !cycleIterator.atEnd();
    cycleIterator.advance()
  ) {
    auto cycleEdges = cycleIterator.getCurrentCycle();

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
      auto a = getAdjacentOfEdge(intersection.front());
      auto b = getAdjacentOfEdge(intersection.back());

      // Find which adjacents are overlapping
      auto indexPair = std::pair<AtomIndexType, AtomIndexType> {
        std::min(a, b),
        std::max(a, b)
      };

      if(
        linksMap.count(indexPair) == 0
        || (
          linksMap.count(indexPair) > 0
          && cycleEdges.size() < links.at(linksMap.at(indexPair)).cycleSequence.size()
        )
      ) {
        LinkInformation newLink;
        newLink.indexPair = indexPair;
        newLink.cycleSequence = centralizeRingIndexSequence(
          makeRingIndexSequence(cycleEdges, graph),
          source
        );

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

  return links;
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

} // namespace MoleculeManip
