#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>

#include "GraphAlgorithms.h"

#include <iostream>

namespace MoleculeManip {

namespace GraphAlgorithms {

using NodeType = BFSVisitors::TreeGenerator::NodeType;

std::shared_ptr<NodeType> makeTree(
  const GraphType& graph,
  const AtomIndexType& source,
  const unsigned& maxDepth
) {
  using ColorMapBase = std::map<
    AtomIndexType,
    boost::default_color_type
  >;

  ColorMapBase colorMap;
  boost::associative_property_map<ColorMapBase> propColorMap(colorMap);
  boost::queue<GraphType::vertex_descriptor> Q;
  auto basePtr = std::make_shared<
    BFSVisitors::TreeGenerator::NodeType
  >(source);

  BFSVisitors::TreeGenerator visitor(
    source,
    maxDepth,
    basePtr
  );

  try {
    boost::breadth_first_visit(
      // The graph to operate on
      graph,
      // The vertex to start with
      source,
      // A queue object to store vertex_descriptors
      Q,
      // The visitor to use
      visitor,
      // A map to store color (state)
      propColorMap
    );
  } catch(BFSVisitors::TreeGenerator::EarlyExit& e) {}

  return basePtr;
}

unsigned numConnectedComponents(const GraphType& graph) {
  std::vector<AtomIndexType> component(boost::num_vertices(graph));

  return boost::connected_components(graph, &component[0]);
}

/* Returns a set of links between source's substituents
 */
std::set<
  std::pair<AtomIndexType, AtomIndexType>
> findSubstituentLinks(
  const GraphType& graph,
  const AtomIndexType& source,
  const std::set<AtomIndexType>& activeSubstituents
) {
  std::set<
    std::pair<AtomIndexType, AtomIndexType>
  > connectedPairs;

  using ColorMapBase = std::map<
    AtomIndexType,
    boost::default_color_type
  >;

  ColorMapBase colorMap;
  boost::associative_property_map<ColorMapBase> propColorMap(colorMap);
  boost::queue<GraphType::vertex_descriptor> Q;

  BFSVisitors::SubstituentLinkSearcher visitor(
    source,
    activeSubstituents,
    connectedPairs // output
  );

  boost::breadth_first_visit(
    // The graph to operate on
    graph,
    // The vertex to start with
    source,
    // A queue object to store vertex_descriptors
    Q,
    // The visitor to use
    visitor,
    // A map to store color (state)
    propColorMap
  );

  return connectedPairs;
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
