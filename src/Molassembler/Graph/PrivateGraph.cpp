/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Graph/PrivateGraph.h"

#include "Molassembler/Molecule/AtomEnvironmentHash.h"
#include "Molassembler/Molecule/MolGraphWriter.h"
#include "Molassembler/Temple/Functional.h"

#include "boost/graph/graphviz.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/graph/graph_utility.hpp"
#include "boost/graph/breadth_first_search.hpp"
#include "boost/graph/biconnected_components.hpp"
#include "boost/graph/connected_components.hpp"
#include "boost/optional.hpp"


namespace Scine {
namespace Molassembler {
namespace {

//! Visitor to help split a graph along a bridge edge
struct BridgeSplittingBFSVisitor {
  using Graph = PrivateGraph::BglType;
  using Vertex = PrivateGraph::Vertex;
  using Edge = PrivateGraph::Edge;
  using BitsetPtr = std::shared_ptr<
    std::vector<bool>
  >;

  BridgeSplittingBFSVisitor() = default;

  BridgeSplittingBFSVisitor(
    AtomIndex left,
    const std::vector<AtomIndex>& site,
    const Graph& graph,
    BitsetPtr ptr
  ) : bitsetPtr(std::move(ptr)) {
    initialized = site;
    initialized.push_back(left);
    *bitsetPtr = std::vector<bool>(boost::num_vertices(graph), false);
    for(const AtomIndex i : site) {
      bitsetPtr->at(i) = true;
    }
  }

  BridgeSplittingBFSVisitor(Edge b, const Graph& graph, BitsetPtr ptr)
    : BridgeSplittingBFSVisitor(
      boost::source(b, graph),
      std::vector<AtomIndex> {1, boost::target(b, graph)},
      graph,
      std::move(ptr)
    ) {}

  void initialize_vertex(Vertex /* v */, const Graph& /* g */) {}
  void discover_vertex(Vertex /* v */, const Graph& /* g */) {}
  void examine_vertex(Vertex /* v */, const Graph& /* g */) {}
  void finish_vertex(Vertex /* v */, const Graph& /* g */) {}
  void examine_edge(const Edge& /* e */, const Graph& /* g */) {}
  void non_tree_edge(const Edge& /* e */, const Graph& /* g */) {}
  void gray_target(const Edge& /* e */, const Graph& /* g */) {}
  void black_target(const Edge& /* e */, const Graph& /* g */) {}

  void tree_edge(const Edge& e, const Graph& g) {
    const Vertex target = boost::target(e, g);
    // Need to make sure not to overwrite the preset values
    if(!Temple::makeContainsPredicate(initialized)(target)) {
      // Copy the source bitset value to the target
      const Vertex source = boost::source(e, g);
      bitsetPtr->at(target) = bitsetPtr->at(source);
    }
  }

  std::vector<AtomIndex> initialized;
  BitsetPtr bitsetPtr;
};

} // namespace

constexpr PrivateGraph::Vertex PrivateGraph::removalPlaceholder;

/* Constructors */
PrivateGraph::PrivateGraph() = default;
PrivateGraph::PrivateGraph(const PrivateGraph::Vertex N) : graph_ {N} {}

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
PrivateGraph::PrivateGraph(const PrivateGraph& other) : graph_(other.graph_) {}
PrivateGraph::PrivateGraph(PrivateGraph&& other) = default;

PrivateGraph& PrivateGraph::operator = (const PrivateGraph& other) {
  graph_ = other.graph_;
  properties_.invalidate();
  return *this;
}
PrivateGraph& PrivateGraph::operator = (PrivateGraph&& other) noexcept {
  graph_ = std::move(other.graph_);
  properties_ = std::move(other.properties_);
  return *this;
}
PrivateGraph::~PrivateGraph() = default;

/* Modifiers */
PrivateGraph::Edge PrivateGraph::addEdge(const Vertex a, const Vertex b, const BondType bondType) {
  /* We have to be careful here since the edge list for a vector is not a set
   * (see BglType).  Check if there is already such an edge before adding it.
   */
  auto existingEdgePair = boost::edge(a, b, graph_);
  if(existingEdgePair.second) {
    throw std::logic_error("Edge already exists!");
  }

  // Invalidate the cache only after all throwing conditions
  properties_.invalidate();

  auto newBondPair = boost::add_edge(a, b, graph_);

  // The bond may not yet exist
  assert(newBondPair.second);
  graph_[newBondPair.first].bondType = bondType;

  return newBondPair.first;
}

PrivateGraph::Vertex PrivateGraph::addVertex(const Utils::ElementType elementType) {
  // Invalidate the cache values
  properties_.invalidate();

  PrivateGraph::Vertex newVertex = boost::add_vertex(graph_);
  graph_[newVertex].elementType = elementType;
  return newVertex;
}

void PrivateGraph::applyPermutation(const std::vector<Vertex>& permutation) {
  // Invalidate the cache
  properties_.invalidate();

  BglType transformedGraph(boost::num_vertices(graph_));

  /* I failed to get copy_graph to perform the permutation for me, so we copy
   * manually:
   */
  for(Vertex i : vertices()) {
    transformedGraph[permutation.at(i)].elementType = graph_[i].elementType;
  }

  for(Edge e : edges()) {
    auto newBondPair = boost::add_edge(
      permutation.at(boost::source(e, graph_)),
      permutation.at(boost::target(e, graph_)),
      transformedGraph
    );

    transformedGraph[newBondPair.first].bondType = graph_[e].bondType;
  }

  std::swap(graph_, transformedGraph);
}

BondType& PrivateGraph::bondType(const PrivateGraph::Edge& edge) {
  // Invalidate the cache values
  properties_.invalidate();

  return graph_[edge].bondType;
}

void PrivateGraph::clearVertex(Vertex a) {
  // Invalidate the cache values
  properties_.invalidate();

  boost::clear_vertex(a, graph_);
}

void PrivateGraph::removeEdge(const Edge& e) {
  // Invalidate the cache values
  properties_.invalidate();

  boost::remove_edge(e, graph_);
}

void PrivateGraph::removeVertex(Vertex a) {
  // Invalidate the cache values
  properties_.invalidate();

  boost::remove_vertex(a, graph_);
}

std::unordered_map<PrivateGraph::Vertex, PrivateGraph::Vertex>
PrivateGraph::merge(
  const PrivateGraph& other,
  const std::vector<Vertex>& copyVertices
) {
  /* Copy over element types, creating the new vertices, collecting
   * the new vertex indices
   */
  std::unordered_map<Vertex, Vertex> copyVertexTargetIndices;
  if(copyVertices.empty()) {
    for(const Vertex vertex : other.vertices()) {
      copyVertexTargetIndices.emplace(
        vertex,
        addVertex(other.elementType(vertex))
      );
    }
  } else {
    for(const Vertex vertexIndex : copyVertices) {
      copyVertexTargetIndices.emplace(
        vertexIndex,
        addVertex(other.elementType(vertexIndex))
      );
    }
  }

  for(const PrivateGraph::Edge& e : other.edges()) {
    const auto findSourceIter = copyVertexTargetIndices.find(
      other.source(e)
    );

    const auto findTargetIter = copyVertexTargetIndices.find(
      other.target(e)
    );

    /* If both source and target indices are keys in copyVertexTargetIndices,
     * then we copy the edge
     */
    if(
      findSourceIter != std::end(copyVertexTargetIndices)
      && findTargetIter != std::end(copyVertexTargetIndices)
    ) {
      addEdge(
        findSourceIter->second,
        findTargetIter->second,
        other.bondType(e)
      );
    }
  }

  return copyVertexTargetIndices;
}

bool PrivateGraph::adjacent(const Vertex a, const Vertex b) const {
  auto edge = boost::edge(a, b, graph_);
  return edge.second;
}

Utils::ElementType& PrivateGraph::elementType(const Vertex a) {
  // Invalidate the cache values
  properties_.invalidate();

  return graph_[a].elementType;
}

PrivateGraph::BglType& PrivateGraph::bgl() {
  // Invalidate the cache values
  properties_.invalidate();

  return graph_;
}

/* Information */

bool PrivateGraph::canRemove(const Vertex a) const {
  assert(a < V());

  /* A molecule is at least one atom. Conceptually, a molecule should consist
   * of at least two atoms, but this is done for usability.
   */
  if(V() == 1) {
    return false;
  }

  // Removable if the vertex is not an articulation vertex
  return removalSafetyData().articulationVertices.count(a) == 0;
}

bool PrivateGraph::canRemove(const Edge& edge) const {
  // Make sure the edge exists in the first place
  assert(
    boost::edge(
      boost::source(edge, graph_),
      boost::target(edge, graph_),
      graph_
    ).second
  );

  const auto& bridges = removalSafetyData().bridges;

  // Removable if the edge is not a bridge
  return bridges.count(edge) == 0;
}

/* Information */
unsigned PrivateGraph::connectedComponents() const {
  std::vector<unsigned> componentMap(V());
  return boost::connected_components(graph_, &componentMap[0]);
}

unsigned PrivateGraph::connectedComponents(std::vector<unsigned>& componentMap) const {
  const Vertex size = V();

  if(componentMap.size() != size) {
    componentMap.resize(size);
  }
  return boost::connected_components(graph_, &componentMap[0]);
}

BondType PrivateGraph::bondType(const PrivateGraph::Edge& edge) const {
  return graph_[edge].bondType;
}

Utils::ElementType PrivateGraph::elementType(const Vertex a) const {
  return graph_[a].elementType;
}

PrivateGraph::Edge PrivateGraph::edge(const Vertex a, const Vertex b) const {
  auto edge = boost::edge(a, b, graph_);

  if(!edge.second) {
    throw std::out_of_range("Specified edge does not exist in the graph");
  }

  return edge.first;
}

boost::optional<PrivateGraph::Edge> PrivateGraph::edgeOption(const Vertex a, const Vertex b) const {
  auto edge = boost::edge(a, b, graph_);

  if(edge.second) {
    return edge.first;
  }

  return boost::none;
}

Utils::ElementTypeCollection PrivateGraph::elementCollection() const {
  const Vertex size = V();

  Utils::ElementTypeCollection elements;
  elements.reserve(size);

  for(AtomIndex i = 0; i < size; ++i) {
    elements.push_back(elementType(i));
  }

  return elements;
}

PrivateGraph::Vertex PrivateGraph::source(const PrivateGraph::Edge& edge) const {
  return boost::source(edge, graph_);
}

PrivateGraph::Vertex PrivateGraph::target(const PrivateGraph::Edge& edge) const {
  return boost::target(edge, graph_);
}

unsigned PrivateGraph::degree(const PrivateGraph::Vertex a) const {
  return boost::out_degree(a, graph_);
}

std::string PrivateGraph::graphviz() const {
  MolGraphWriter propertyWriter(&*this, nullptr);

  std::stringstream graphvizStream;

  boost::write_graphviz(
    graphvizStream,
    bgl(),
    propertyWriter,
    propertyWriter,
    propertyWriter
  );

  return graphvizStream.str();
}

PrivateGraph::Vertex PrivateGraph::N() const {
  return boost::num_vertices(graph_);
}

PrivateGraph::Vertex PrivateGraph::B() const {
  return boost::num_edges(graph_);
}

PrivateGraph::Vertex PrivateGraph::V() const {
  return boost::num_vertices(graph_);
}

unsigned PrivateGraph::E() const {
  return boost::num_edges(graph_);
}

boost::optional<std::vector<AtomIndex>> PrivateGraph::modularIsomorphism(
  const PrivateGraph& other,
  const AtomEnvironmentComponents components
) const {
  const unsigned thisNumAtoms = V();

  //! Quick checks
  if(thisNumAtoms != other.V() || E() != other.E()) {
    return boost::none;
  }

  std::vector<Hashes::HashType> thisHashes;
  std::vector<Hashes::HashType> otherHashes;
  Hashes::HashType maxHash;
  std::tie(thisHashes, otherHashes, maxHash) = Hashes::narrow(
    Hashes::generate(*this, boost::none, components),
    Hashes::generate(other, boost::none, components)
  );

  std::vector<AtomIndex> indexMap(thisNumAtoms);

  const bool isomorphic = boost::isomorphism(
    bgl(),
    other.bgl(),
    boost::make_safe_iterator_property_map(
      indexMap.begin(),
      thisNumAtoms,
      boost::get(boost::vertex_index, bgl())
    ),
    Hashes::LookupFunctor(thisHashes),
    Hashes::LookupFunctor(otherHashes),
    maxHash,
    boost::get(boost::vertex_index, bgl()),
    boost::get(boost::vertex_index, other.bgl())
  );

  if(isomorphic) {
    return indexMap;
  }

  return boost::none;
}

bool PrivateGraph::identicalGraph(const PrivateGraph& other) const {
  assert(V() == other.V() && E() == other.E());

  // Make sure topology matches
  return Temple::all_of(
    edges(),
    [&](const Edge& edge) -> bool {
      bool edgeExists;

      std::tie(std::ignore, edgeExists) = boost::edge(
        boost::source(edge, graph_),
        boost::target(edge, graph_),
        other.graph_
      );

      return edgeExists;
    }
  );
}

std::pair<
  std::vector<AtomIndex>,
  std::vector<AtomIndex>
> PrivateGraph::splitAlongBridge(Edge bridge) const {
  return splitAlongBridge(
    source(bridge),
    std::vector<AtomIndex>(1, target(bridge))
  );
}

std::pair<
  std::vector<AtomIndex>,
  std::vector<AtomIndex>
> PrivateGraph::splitAlongBridge(AtomIndex left, const std::vector<AtomIndex>& right) const {
  // Can't really check that it's truly a bridge edge since there's a bunch of edges
  if(right.empty()) {
    throw std::invalid_argument("The atom set on the right may not be empty");
  }

  if(right.size() == 1) {
    const auto maybeEdge = edgeOption(left, right.front());
    if(!maybeEdge) {
      throw std::invalid_argument("Left and right atom are not connected");
    }
    if(removalSafetyData().bridges.count(maybeEdge.value()) == 0) {
      throw std::invalid_argument("The supplied edge is not a bridge edge");
    }
  } else {
    /* Can't really check if the edges are truly collectively a bridge edge
     * without making a copy of the graph to check
     */
    if(
      Temple::any_of(
        right,
        [&](const AtomIndex r) -> bool {
          const auto maybeEdge = edgeOption(left, r);
          return !maybeEdge || bondType(maybeEdge.value()) != BondType::Eta;
        }
      )
    ) {
      throw std::invalid_argument("Right atom set is not a haptic site of the left atom");
    }
  }

  auto bitsetPtr = std::make_shared<std::vector<bool>>();
  BridgeSplittingBFSVisitor visitor(left, right, graph_, bitsetPtr);

  boost::breadth_first_search(
    graph_,
    right.front(),
    boost::visitor(visitor)
  );

  // Transform the bitset into a left and right
  std::vector<AtomIndex> lefts;
  std::vector<AtomIndex> rights;
  const unsigned size = V();
  lefts.reserve(size);
  rights.reserve(size);

  for(AtomIndex i = 0; i < size; ++i) {
    if(bitsetPtr->at(i)) {
      rights.push_back(i);
    } else {
      lefts.push_back(i);
    }
  }

  lefts.shrink_to_fit();
  rights.shrink_to_fit();

  return std::make_pair(
    std::move(lefts),
    std::move(rights)
  );
}

PrivateGraph::VertexRange PrivateGraph::vertices() const {
  auto iters = boost::vertices(graph_);
  return {
    std::move(iters.first),
    std::move(iters.second)
  };
}

PrivateGraph::EdgeRange PrivateGraph::edges() const {
  auto iters = boost::edges(graph_);
  return {
    std::move(iters.first),
    std::move(iters.second)
  };
}

PrivateGraph::AdjacentVertexRange PrivateGraph::adjacents(const Vertex a) const {
  auto iters = boost::adjacent_vertices(a, graph_);
  return {
    std::move(iters.first),
    std::move(iters.second)
  };
}

PrivateGraph::IncidentEdgeRange PrivateGraph::edges(const Vertex a) const {
  auto iters = boost::out_edges(a, graph_);
  return {
    std::move(iters.first),
    std::move(iters.second)
  };
}

bool PrivateGraph::operator == (const PrivateGraph& other) const {
  // Better with newer boost: has_value()
  return static_cast<bool>(
    modularIsomorphism(other, AtomEnvironmentComponents::All)
  );
}

const PrivateGraph::BglType& PrivateGraph::bgl() const {
  return graph_;
}

void PrivateGraph::populateProperties() const {
  if(!properties_.removalSafetyDataOption) {
    properties_.removalSafetyDataOption = generateRemovalSafetyData_();
  }
  if(!properties_.cyclesOption) {
    properties_.cyclesOption = generateCycles_();
  }
}

const PrivateGraph::RemovalSafetyData& PrivateGraph::removalSafetyData() const {
  if(!properties_.removalSafetyDataOption) {
    properties_.removalSafetyDataOption = generateRemovalSafetyData_();
  }

  return *properties_.removalSafetyDataOption;
}

const Cycles& PrivateGraph::cycles() const {
  if(!properties_.cyclesOption) {
    properties_.cyclesOption = generateCycles_();
  }

  return *properties_.cyclesOption;
}

const Cycles& PrivateGraph::etaPreservedCycles() const {
  if(!properties_.etaPreservedCyclesOption) {
    properties_.etaPreservedCyclesOption = generateEtaPreservedCycles_();
  }

  return *properties_.etaPreservedCyclesOption;
}

PrivateGraph::RemovalSafetyData PrivateGraph::generateRemovalSafetyData_() const {
  RemovalSafetyData safetyData;

  std::vector<PrivateGraph::Vertex> articulationVertices;

  using ComponentMapBase = std::map<PrivateGraph::Edge, std::size_t>;

  ComponentMapBase componentMapData;
  boost::associative_property_map<ComponentMapBase> componentMap(componentMapData);
  std::size_t numComponents;

  // Calculate the biconnected components and articulation vertices
  std::tie(numComponents, std::ignore) = boost::biconnected_components(
    graph_,
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

Cycles PrivateGraph::generateCycles_() const {
  return {*this};
}

Cycles PrivateGraph::generateEtaPreservedCycles_() const {
  return {*this, false};
}

} // namespace Molassembler
} // namespace Scine
