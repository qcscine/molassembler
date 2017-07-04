#include "GraphAlgorithms.h"

#include <boost/graph/connected_components.hpp>

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

std::shared_ptr<NodeType> makeTree(
  const GraphType& graph,
  const AtomIndexType& source
) {
  return makeTree(
    graph,
    source,
    0
  );
}

std::shared_ptr<NodeType> makeTree(
  const GraphType& graph
) {
  return makeTree(
    graph,
    0,
    0
  );
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

} // namespace GraphAlgorithms

} // namespace MoleculeManip
