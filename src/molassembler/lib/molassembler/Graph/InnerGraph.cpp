/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Graph/InnerGraph.h"

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
  using Graph = InnerGraph::BGLType;
  using Vertex = InnerGraph::Vertex;
  using Edge = InnerGraph::Edge;

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

constexpr InnerGraph::Vertex InnerGraph::removalPlaceholder;

/* Constructors */
InnerGraph::InnerGraph() = default;
InnerGraph::InnerGraph(const InnerGraph::Vertex N) : _graph {N} {}

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
InnerGraph::InnerGraph(const InnerGraph& other) : _graph(other._graph) {}
InnerGraph::InnerGraph(InnerGraph&& other) = default;
InnerGraph& InnerGraph::operator = (const InnerGraph& other) {
  _graph = other._graph;
  _properties.invalidate();
  return *this;
}
InnerGraph& InnerGraph::operator = (InnerGraph&& other) {
  _graph = std::move(other._graph);
  _properties = std::move(other._properties);
  return *this;
}
InnerGraph::~InnerGraph() = default;

/* Modifiers */
InnerGraph::Edge InnerGraph::addEdge(const Vertex a, const Vertex b, const BondType bondType) {
  /* We have to be careful here since the edge list for a vector is not a set
   * (see BGLType).  Check if there is already such an edge before adding it.
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

InnerGraph::Vertex InnerGraph::addVertex(const Utils::ElementType elementType) {
  // Invalidate the cache values
  _properties.invalidate();

  InnerGraph::Vertex newVertex = boost::add_vertex(_graph);
  _graph[newVertex].elementType = elementType;
  return newVertex;
}

void InnerGraph::applyPermutation(const std::vector<Vertex>& permutation) {
  // Invalidate the cache
  _properties.invalidate();

  BGLType transformedGraph(boost::num_vertices(_graph));

  /* I failed to get copy_graph to perform the permutation for me, so we copy
   * manually:
   */
  for(Vertex i : boost::make_iterator_range(vertices())) {
    transformedGraph[permutation.at(i)].elementType = _graph[i].elementType;
  }

  for(Edge e : boost::make_iterator_range(edges())) {
    auto newBondPair = boost::add_edge(
      permutation.at(boost::source(e, _graph)),
      permutation.at(boost::target(e, _graph)),
      transformedGraph
    );

    transformedGraph[newBondPair.first].bondType = _graph[e].bondType;
  }

  std::swap(_graph, transformedGraph);
}

BondType& InnerGraph::bondType(const InnerGraph::Edge& edge) {
  // Invalidate the cache values
  _properties.invalidate();

  return _graph[edge].bondType;
}

void InnerGraph::clearVertex(Vertex a) {
  // Invalidate the cache values
  _properties.invalidate();

  boost::clear_vertex(a, _graph);
}

void InnerGraph::removeEdge(const Edge& e) {
  // Invalidate the cache values
  _properties.invalidate();

  boost::remove_edge(e, _graph);
}

void InnerGraph::removeVertex(Vertex a) {
  // Invalidate the cache values
  _properties.invalidate();

  boost::remove_vertex(a, _graph);
}

Utils::ElementType& InnerGraph::elementType(const Vertex a) {
  // Invalidate the cache values
  _properties.invalidate();

  return _graph[a].elementType;
}

InnerGraph::BGLType& InnerGraph::bgl() {
  // Invalidate the cache values
  _properties.invalidate();

  return _graph;
}

/* Information */

bool InnerGraph::canRemove(const Vertex a) const {
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

bool InnerGraph::canRemove(const Edge& edge) const {
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
unsigned InnerGraph::connectedComponents() const {
  std::vector<unsigned> componentMap(N());
  return boost::connected_components(_graph, &componentMap[0]);
}

unsigned InnerGraph::connectedComponents(std::vector<unsigned>& componentMap) const {
  const Vertex size = N();

  if(componentMap.size() != size) {
    componentMap.resize(size);
  }
  return boost::connected_components(_graph, &componentMap[0]);
}

BondType InnerGraph::bondType(const InnerGraph::Edge& edge) const {
  return _graph[edge].bondType;
}

Utils::ElementType InnerGraph::elementType(const Vertex a) const {
  return _graph[a].elementType;
}

InnerGraph::Edge InnerGraph::edge(const Vertex a, const Vertex b) const {
  auto edge = boost::edge(a, b, _graph);
  assert(edge.second);
  return edge.first;
}

boost::optional<InnerGraph::Edge> InnerGraph::edgeOption(const Vertex a, const Vertex b) const {
  auto edge = boost::edge(a, b, _graph);

  if(edge.second) {
    return edge.first;
  }

  return boost::none;
}

InnerGraph::Vertex InnerGraph::source(const InnerGraph::Edge& edge) const {
  return boost::source(edge, _graph);
}

InnerGraph::Vertex InnerGraph::target(const InnerGraph::Edge& edge) const {
  return boost::target(edge, _graph);
}

InnerGraph::Vertex InnerGraph::degree(const InnerGraph::Vertex a) const {
  return boost::out_degree(a, _graph);
}

InnerGraph::Vertex InnerGraph::N() const {
  return boost::num_vertices(_graph);
}

InnerGraph::Vertex InnerGraph::B() const {
  return boost::num_edges(_graph);
}

bool InnerGraph::identicalGraph(const InnerGraph& other) const {
  assert(N() == other.N() && B() == other.B());

  // Make sure topology matches
  return temple::all_of(
    boost::make_iterator_range(edges()),
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
> InnerGraph::splitAlongBridge(Edge bridge) const {
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

InnerGraph::VertexRange InnerGraph::vertices() const {
  return boost::vertices(_graph);
}

InnerGraph::EdgeRange InnerGraph::edges() const {
  return boost::edges(_graph);
}

InnerGraph::AdjacentVertexRange InnerGraph::adjacents(const Vertex a) const {
  return boost::adjacent_vertices(a, _graph);
}

InnerGraph::IncidentEdgeRange InnerGraph::edges(const Vertex a) const {
  return boost::out_edges(a, _graph);
}

const InnerGraph::BGLType& InnerGraph::bgl() const {
  return _graph;
}

void InnerGraph::populateProperties() const {
  if(!_properties.removalSafetyDataOption) {
    _properties.removalSafetyDataOption = _generateRemovalSafetyData();
  }
  if(!_properties.cyclesOption) {
    _properties.cyclesOption = _generateCycles();
  }
}

const InnerGraph::RemovalSafetyData& InnerGraph::removalSafetyData() const {
  if(!_properties.removalSafetyDataOption) {
    _properties.removalSafetyDataOption = _generateRemovalSafetyData();
  }

  return *_properties.removalSafetyDataOption;
}

const Cycles& InnerGraph::cycles() const {
  if(!_properties.cyclesOption) {
    _properties.cyclesOption = _generateCycles();
  }

  return *_properties.cyclesOption;
}

const Cycles& InnerGraph::etaPreservedCycles() const {
  if(!_properties.etaPreservedCyclesOption) {
    _properties.etaPreservedCyclesOption = _generateEtaPreservedCycles();
  }

  return *_properties.etaPreservedCyclesOption;
}

InnerGraph::RemovalSafetyData InnerGraph::_generateRemovalSafetyData() const {
  RemovalSafetyData safetyData;

  std::vector<InnerGraph::Vertex> articulationVertices;

  using ComponentMapBase = std::map<InnerGraph::Edge, std::size_t>;

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
    std::set<InnerGraph::Edge>
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

Cycles InnerGraph::_generateCycles() const {
  return Cycles(*this);
}

Cycles InnerGraph::_generateEtaPreservedCycles() const {
  return Cycles(*this, false);
}

} // namespace molassembler

} // namespace Scine
