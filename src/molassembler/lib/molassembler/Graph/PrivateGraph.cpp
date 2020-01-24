/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Graph/PrivateGraph.h"

#include "boost/graph/breadth_first_search.hpp"
#include "boost/graph/biconnected_components.hpp"
#include "boost/graph/connected_components.hpp"
#include "boost/optional.hpp"

#include "temple/Functional.h"

namespace Scine {

namespace molassembler {

namespace detail {

//! Visitor to help split a graph along a bridge edge
struct BridgeSplittingBFSVisitor {
  using Graph = PrivateGraph::BglType;
  using Vertex = PrivateGraph::Vertex;
  using Edge = PrivateGraph::Edge;

  AtomIndex left, right;
  using BitsetPtr = std::shared_ptr<
    std::vector<bool>
  >;

  BitsetPtr bitsetPtr;

  BridgeSplittingBFSVisitor() = default;
  BridgeSplittingBFSVisitor(Edge b, const Graph& graph, BitsetPtr ptr)
    : left(boost::source(b, graph)),
      right(boost::target(b, graph)),
      bitsetPtr(std::move(ptr))
  {
    *bitsetPtr = std::vector<bool>(boost::num_vertices(graph), false);
    bitsetPtr->at(right) = true;
  }

  void initialize_vertex(Vertex /* v */, const Graph& /* g */) {}
  void discover_vertex(Vertex /* v */, const Graph& /* g */) {}
  void examine_vertex(Vertex /* v */, const Graph& /* g */) {}
  void finish_vertex(Vertex /* v */, const Graph& /* g */) {}
  void examine_edge(const Edge& /* e */, const Graph& /* g */) {}

  void tree_edge(const Edge& e, const Graph& g) {
    Vertex target = boost::target(e, g);
    if(target != left) {
      // Copy the source bitset value to the target
      bitsetPtr->at(target) = bitsetPtr->at(
        boost::source(e, g)
      );
    }
  }

  void non_tree_edge(const Edge& /* e */, const Graph& /* g */) {}
  void gray_target(const Edge& /* e */, const Graph& /* g */) {}
  void black_target(const Edge& /* e */, const Graph& /* g */) {}
};

} // namespace detail

constexpr PrivateGraph::Vertex PrivateGraph::removalPlaceholder;

/* Constructors */
PrivateGraph::PrivateGraph() = default;
PrivateGraph::PrivateGraph(const PrivateGraph::Vertex N) : _graph {N} {}

/* Rule of five members */
/* Implementation note
 *
 * Copying the graph invalidates any stored edge_descriptors in the properties
 * since the copied graph is now elsewhere in memory and edge_descriptor's
 * property member is no longer correct, leading to subtle bugs.
 *
 * Moving the graph does not invalidate descriptors since the graph does
 * not change location in memory.
 */
PrivateGraph::PrivateGraph(const PrivateGraph& other) : _graph(other._graph) {}
PrivateGraph::PrivateGraph(PrivateGraph&& other) = default;
PrivateGraph& PrivateGraph::operator = (const PrivateGraph& other) {
  _graph = other._graph;
  _properties.invalidate();
  return *this;
}
PrivateGraph& PrivateGraph::operator = (PrivateGraph&& other) {
  _graph = std::move(other._graph);
  _properties = std::move(other._properties);
  return *this;
}
PrivateGraph::~PrivateGraph() = default;

/* Modifiers */
PrivateGraph::Edge PrivateGraph::addEdge(const Vertex a, const Vertex b, const BondType bondType) {
  /* We have to be careful here since the edge list for a vector is not a set
   * (see BglType).  Check if there is already such an edge before adding it.
   */
  auto existingEdgePair = boost::edge(a, b, _graph);
  if(existingEdgePair.second) {
    throw std::logic_error("Edge already exists!");
  }

  // Invalidate the cache only after all throwing conditions
  _properties.invalidate();

  auto newBondPair = boost::add_edge(a, b, _graph);

  // The bond may not yet exist
  assert(newBondPair.second);
  _graph[newBondPair.first].bondType = bondType;

  return newBondPair.first;
}

PrivateGraph::Vertex PrivateGraph::addVertex(const Utils::ElementType elementType) {
  // Invalidate the cache values
  _properties.invalidate();

  PrivateGraph::Vertex newVertex = boost::add_vertex(_graph);
  _graph[newVertex].elementType = elementType;
  return newVertex;
}

void PrivateGraph::applyPermutation(const std::vector<Vertex>& permutation) {
  // Invalidate the cache
  _properties.invalidate();

  BglType transformedGraph(boost::num_vertices(_graph));

  /* I failed to get copy_graph to perform the permutation for me, so we copy
   * manually:
   */
  for(Vertex i : vertices()) {
    transformedGraph[permutation.at(i)].elementType = _graph[i].elementType;
  }

  for(Edge e : edges()) {
    auto newBondPair = boost::add_edge(
      permutation.at(boost::source(e, _graph)),
      permutation.at(boost::target(e, _graph)),
      transformedGraph
    );

    transformedGraph[newBondPair.first].bondType = _graph[e].bondType;
  }

  std::swap(_graph, transformedGraph);
}

BondType& PrivateGraph::bondType(const PrivateGraph::Edge& edge) {
  // Invalidate the cache values
  _properties.invalidate();

  return _graph[edge].bondType;
}

void PrivateGraph::clearVertex(Vertex a) {
  // Invalidate the cache values
  _properties.invalidate();

  boost::clear_vertex(a, _graph);
}

void PrivateGraph::removeEdge(const Edge& e) {
  // Invalidate the cache values
  _properties.invalidate();

  boost::remove_edge(e, _graph);
}

void PrivateGraph::removeVertex(Vertex a) {
  // Invalidate the cache values
  _properties.invalidate();

  boost::remove_vertex(a, _graph);
}

Utils::ElementType& PrivateGraph::elementType(const Vertex a) {
  // Invalidate the cache values
  _properties.invalidate();

  return _graph[a].elementType;
}

PrivateGraph::BglType& PrivateGraph::bgl() {
  // Invalidate the cache values
  _properties.invalidate();

  return _graph;
}

/* Information */

bool PrivateGraph::canRemove(const Vertex a) const {
  assert(a < N());

  /* A molecule is at least one atom. Conceptually, a molecule should consist
   * of at least two atoms, but this is done for usability.
   */
  if(N() == 1) {
    return false;
  }

  // Removable if the vertex is not an articulation vertex
  return removalSafetyData().articulationVertices.count(a) == 0;
}

bool PrivateGraph::canRemove(const Edge& edge) const {
  // Make sure the edge exists in the first place
  assert(
    boost::edge(
      boost::source(edge, _graph),
      boost::target(edge, _graph),
      _graph
    ).second
  );

  const auto& bridges = removalSafetyData().bridges;

  // Removable if the edge is not a bridge
  return bridges.count(edge) == 0;
}

/* Information */
unsigned PrivateGraph::connectedComponents() const {
  std::vector<unsigned> componentMap(N());
  return boost::connected_components(_graph, &componentMap[0]);
}

unsigned PrivateGraph::connectedComponents(std::vector<unsigned>& componentMap) const {
  const Vertex size = N();

  if(componentMap.size() != size) {
    componentMap.resize(size);
  }
  return boost::connected_components(_graph, &componentMap[0]);
}

BondType PrivateGraph::bondType(const PrivateGraph::Edge& edge) const {
  return _graph[edge].bondType;
}

Utils::ElementType PrivateGraph::elementType(const Vertex a) const {
  return _graph[a].elementType;
}

PrivateGraph::Edge PrivateGraph::edge(const Vertex a, const Vertex b) const {
  auto edge = boost::edge(a, b, _graph);

  if(!edge.second) {
    throw std::out_of_range("Specified edge does not exist in the graph");
  }

  return edge.first;
}

boost::optional<PrivateGraph::Edge> PrivateGraph::edgeOption(const Vertex a, const Vertex b) const {
  auto edge = boost::edge(a, b, _graph);

  if(edge.second) {
    return edge.first;
  }

  return boost::none;
}

PrivateGraph::Vertex PrivateGraph::source(const PrivateGraph::Edge& edge) const {
  return boost::source(edge, _graph);
}

PrivateGraph::Vertex PrivateGraph::target(const PrivateGraph::Edge& edge) const {
  return boost::target(edge, _graph);
}

PrivateGraph::Vertex PrivateGraph::degree(const PrivateGraph::Vertex a) const {
  return boost::out_degree(a, _graph);
}

PrivateGraph::Vertex PrivateGraph::N() const {
  return boost::num_vertices(_graph);
}

PrivateGraph::Vertex PrivateGraph::B() const {
  return boost::num_edges(_graph);
}

bool PrivateGraph::identicalGraph(const PrivateGraph& other) const {
  assert(N() == other.N() && B() == other.B());

  // Make sure topology matches
  return temple::all_of(
    edges(),
    [&](const Edge& edge) -> bool {
      Edge correspondingEdge;
      bool edgeExists;
      std::tie(correspondingEdge, edgeExists) = boost::edge(
        boost::source(edge, _graph),
        boost::target(edge, _graph),
        other._graph
      );

      return edgeExists;
    }
  );
}

std::pair<
  std::vector<AtomIndex>,
  std::vector<AtomIndex>
> PrivateGraph::splitAlongBridge(Edge bridge) const {
  if(removalSafetyData().bridges.count(bridge) == 0) {
    throw std::invalid_argument("The supplied edge is not a bridge edge");
  }

  auto bitsetPtr = std::make_shared<
    std::vector<bool>
  >();
  detail::BridgeSplittingBFSVisitor visitor(bridge, _graph, bitsetPtr);

  boost::breadth_first_search(
    _graph,
    target(bridge),
    boost::visitor(visitor)
  );

  // Transform the bitset into a left and right
  std::vector<AtomIndex> left, right;
  left.reserve(N());
  right.reserve(N());

  for(AtomIndex i = 0; i < N(); ++i) {
    if(bitsetPtr->at(i)) {
      right.push_back(i);
    } else {
      left.push_back(i);
    }
  }

  left.shrink_to_fit();
  right.shrink_to_fit();

  return std::make_pair(
    std::move(left),
    std::move(right)
  );
}

PrivateGraph::VertexRange PrivateGraph::vertices() const {
  auto iters = boost::vertices(_graph);
  return {
    std::move(iters.first),
    std::move(iters.second)
  };
}

PrivateGraph::EdgeRange PrivateGraph::edges() const {
  auto iters = boost::edges(_graph);
  return {
    std::move(iters.first),
    std::move(iters.second)
  };
}

PrivateGraph::AdjacentVertexRange PrivateGraph::adjacents(const Vertex a) const {
  auto iters = boost::adjacent_vertices(a, _graph);
  return {
    std::move(iters.first),
    std::move(iters.second)
  };
}

PrivateGraph::IncidentEdgeRange PrivateGraph::edges(const Vertex a) const {
  auto iters = boost::out_edges(a, _graph);
  return {
    std::move(iters.first),
    std::move(iters.second)
  };
}

const PrivateGraph::BglType& PrivateGraph::bgl() const {
  return _graph;
}

void PrivateGraph::populateProperties() const {
  if(!_properties.removalSafetyDataOption) {
    _properties.removalSafetyDataOption = _generateRemovalSafetyData();
  }
  if(!_properties.cyclesOption) {
    _properties.cyclesOption = _generateCycles();
  }
}

const PrivateGraph::RemovalSafetyData& PrivateGraph::removalSafetyData() const {
  if(!_properties.removalSafetyDataOption) {
    _properties.removalSafetyDataOption = _generateRemovalSafetyData();
  }

  return *_properties.removalSafetyDataOption;
}

const Cycles& PrivateGraph::cycles() const {
  if(!_properties.cyclesOption) {
    _properties.cyclesOption = _generateCycles();
  }

  return *_properties.cyclesOption;
}

const Cycles& PrivateGraph::etaPreservedCycles() const {
  if(!_properties.etaPreservedCyclesOption) {
    _properties.etaPreservedCyclesOption = _generateEtaPreservedCycles();
  }

  return *_properties.etaPreservedCyclesOption;
}

PrivateGraph::RemovalSafetyData PrivateGraph::_generateRemovalSafetyData() const {
  RemovalSafetyData safetyData;

  std::vector<PrivateGraph::Vertex> articulationVertices;

  using ComponentMapBase = std::map<PrivateGraph::Edge, std::size_t>;

  ComponentMapBase componentMapData;
  boost::associative_property_map<ComponentMapBase> componentMap(componentMapData);
  std::size_t numComponents;

  // Calculate the biconnected components and articulation vertices
  std::tie(numComponents, std::ignore) = boost::biconnected_components(
    _graph,
    componentMap,
    std::back_inserter(articulationVertices)
  );

  // Copy articulation vertices to the set
  for(const auto& vertex : articulationVertices) {
    safetyData.articulationVertices.insert(vertex);
  }

  /* Work out from the biconnected components which edges are bridges: If a
   * biconnected component contains only a single edge, that edge is a bridge
   */
  std::vector<
    std::set<PrivateGraph::Edge>
  > componentSets (numComponents);

  for(const auto& mapIterPair : componentMapData) {
    const auto& edgeIndex = mapIterPair.first;
    const auto& componentIndex = mapIterPair.second;

    componentSets.at(componentIndex).insert(edgeIndex);
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

Cycles PrivateGraph::_generateCycles() const {
  return Cycles(*this);
}

Cycles PrivateGraph::_generateEtaPreservedCycles() const {
  return Cycles(*this, false);
}

} // namespace molassembler

} // namespace Scine
