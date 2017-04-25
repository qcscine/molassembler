#include "AdjacencyListAlgorithms.h"

#include <boost/graph/connected_components.hpp>

#include <iostream>

namespace MoleculeManip {

namespace AdjacencyListAlgorithms {

using NodeType = BFSVisitors::TreeGenerator::NodeType;

std::shared_ptr<NodeType> makeTree(
  const AdjacencyList& adjacencies,
  const AtomIndexType& startingFrom,
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
  >(startingFrom);

  BFSVisitors::TreeGenerator visitor(
    startingFrom,
    maxDepth,
    basePtr
  );

  try {
    boost::breadth_first_visit(
      // The graph to operate on
      adjacencies.access(),  
      // The vertex to start with
      startingFrom,
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

std::shared_ptr<NodeType> makeTree(
  const AdjacencyList& adjacencies,
  const AtomIndexType& startingFrom
) {
  return makeTree(
    adjacencies,
    startingFrom,
    0
  );
}

std::shared_ptr<NodeType> makeTree(
  const AdjacencyList& adjacencies
) {
  return makeTree(
    adjacencies,
    0,
    0
  );
}

unsigned numConnectedComponents(const AdjacencyList& adjacencies) {
  auto& graph = adjacencies.access();
  std::vector<AtomIndexType> component(boost::num_vertices(graph));

  return boost::connected_components(graph, &component[0]);
}

} // namespace AdjacencyListAlgorithms

} // namespace MoleculeManip
