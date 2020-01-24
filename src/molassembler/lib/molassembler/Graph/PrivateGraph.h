/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Boost Graph Library wrapper to help in concealing underlying type
 */

#ifndef INCLUDE_MOLASSEMBLER_INNER_GRAPH_H
#define INCLUDE_MOLASSEMBLER_INNER_GRAPH_H

#include "boost/graph/adjacency_list.hpp"
#include "Utils/Geometry/ElementTypes.h"

#include "molassembler/Cycles.h"

#include <limits>

namespace Scine {
namespace molassembler {

/**
 * @brief Library internal graph class wrapping BGL types
 */
class PrivateGraph {
public:
//!@name Member types
//!@{
  //! Information stored at each graph vertex
  struct VertexData {
    Utils::ElementType elementType;
  };

  //! Information stored at each graph edge
  struct EdgeData {
    BondType bondType;
  };

  /*! The type of the BGL graph.
   *
   * An adjacency list is used due to sparsity of the molecular graph.
   */
  using BglType = boost::adjacency_list<
    /* OutEdgeListS = Type of Container for edges of a vertex
     * Options: vector, list, slist, set, multiset, unordered_set
     * Choice: vecS, we have to rigorously test that no parallel edges are
     *   created
     */
    boost::vecS,
    /* VertexListS = Type of Container for vertices
     * Options: vector, list, slist, set, multiset, unordered_set
     * Choice: vecS, removing vertices is rare, keep memory use limited
     * Consequence: operation remove_vertex() invalidates:
     *   - Vertex descriptors / iterators
     *   - Edge descriptors / iterators
     *   - Adjacency iterators
     */
    boost::vecS,
    // Molecular graphs are undirected
    boost::undirectedS,
    // VertexProperty: What information is stored about vertices?
    VertexData,
    // EdgeProperty: What information is stored about edges?
    EdgeData
    // GraphProperty: Omitted, defaults
    // EdgeListS: Omitted, defaults
  >;

  using Vertex = BglType::vertex_descriptor;
  using Edge = BglType::edge_descriptor;

  using VertexRange = IteratorRange<BglType::vertex_iterator>;
  using EdgeRange = IteratorRange<BglType::edge_iterator>;
  using AdjacentVertexRange = IteratorRange<BglType::adjacency_iterator>;
  using IncidentEdgeRange = IteratorRange<BglType::out_edge_iterator>;

  /*!
   * @brief Data class to return removal safety information on the graph
   *
   * @note These are grouped because it is fairly easy to work out bridges from
   * articulation vertices and hence their joint calculation makes a lot of
   * sense. It might also be easier to cache this way without the dependence.
   */
  struct RemovalSafetyData {
    //! Articulation vertices cannot be removed without disconnecting the graph
    std::unordered_set<PrivateGraph::Vertex> articulationVertices;
    //! Bridges are edges that cannot be removed without disconnecting the graph
    std::set<PrivateGraph::Edge> bridges;
  };
//!@}

//!@name Constructors
//!@{
  //! Empty constructor
  PrivateGraph();
  //! Preallocating constructor
  PrivateGraph(Vertex N);
//!@}

//!@name Rule of five members
//!@{
  PrivateGraph(const PrivateGraph& other);
  PrivateGraph(PrivateGraph&& other);
  PrivateGraph& operator = (const PrivateGraph& other);
  PrivateGraph& operator = (PrivateGraph&& other);
  ~PrivateGraph();
//!@}

//!@name Static members
//!@{
  //! Placeholder vertex index used to indicate that a vertex has been removed
  static constexpr Vertex removalPlaceholder = std::numeric_limits<Vertex>::max();
//!@}

//!@name Modification
//!@{
  /*! @brief Add an edge to the graph
   *
   * @complexity{@math{\Theta(1)}}
   * @throws std::logic_error if the edge already exists
   */
  Edge addEdge(Vertex a, Vertex b, BondType bondType);

  /*! @brief Add a (disconnected) vertex to the graph
   *
   * @complexity{@math{\Theta(1)}}
   */
  Vertex addVertex(Utils::ElementType elementType);

  /*! @brief Apply a permutation to the graph
   *
   * @complexity{@math{\Theta(V)}}
   */
  void applyPermutation(const std::vector<Vertex>& permutation);

  /** @brief Fetch the bond type of an edge
   *
   * @complexity{@math{\Theta(1)}}
   * @pre The edge exists
   */
  BondType& bondType(const Edge& edge);

  //! @brief Referential access to the underlying BGL graph
  BglType& bgl();

  /*! @brief Removes all bonds involving a vertex
   *
   * @complexity{@math{\Theta(1)}}
   */
  void clearVertex(Vertex a);

  /*! @brief Fetches the element type of a vertex
   *
   * @complexity{@math{\Theta(1)}}
   * @pre The vertex exists
   */
  Utils::ElementType& elementType(Vertex a);

  /*! @brief Removes an edge from the graph.
   *
   * @complexity{@math{\Theta(1)}}
   */
  void removeEdge(const Edge& e);

  /*! @brief Removes a vertex from the graph.
   *
   * @complexity{@math{\Theta(1)}}
   * @warning This invalidates ALL vertex and edge descriptors!
   */
  void removeVertex(Vertex a);
//!@}

//!@name Information
//!@{
  /*! @brief Fetches the element type of a vertex
   *
   * @complexity{@math{\Theta(1)}}
   * @pre The vertex exists
   */
  Utils::ElementType elementType(Vertex a) const;
  /** @brief Fetch the bond type of an edge
   *
   * @complexity{@math{\Theta(1)}}
   * @pre The edge exists
   */
  BondType bondType(const Edge& edge) const;
  //! @brief Nonmodifiable access to the underlying BGL graph
  const BglType& bgl() const;

  /*! @brief Determine whether a vertex can be safely removed
   *
   * A Vertex can be safely removed if it is not an articulation vertex or if
   * the number of vertices is more than one.
   *
   * @complexity{@math{O(V)} worst case, if removal data is cached
   * @math{\Theta(1)}}
   *
   * @note This function is not thread-safe.
   */
  bool canRemove(Vertex a) const;
  /*! @brief Determine whether an edge can be safely removed
   * An edge can be safely removed if it is not a bridge edge
   *
   * @complexity{@math{O(V)} worst case, if removal data is cached
   * @math{\Theta(1)}}
   *
   * @note This function is not thread-safe.
   */
  bool canRemove(const Edge& edge) const;

  /*! @brief Number of connected components
   *
   * @complexity{@math{\Theta(V)}}
   */
  unsigned connectedComponents() const;

  /*! @brief Connected components, yielding a map from vertices to component index
   *
   * @complexity{@math{\Theta(V)}}
   */
  unsigned connectedComponents(std::vector<unsigned>& componentMap) const;

  /*! @brief Make an edge descriptor from two vertex descriptors
   *
   * @throws std::out_of_range if the edge doesn't exist in the graph
   * @complexity{@math{\Theta(1)}}
   */
  Edge edge(Vertex a, Vertex b) const;

  /*! @brief Get an edge descriptor from two vertex descriptors, get None if the edge doesn't exist
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<Edge> edgeOption(Vertex a, Vertex b) const;

  Vertex source(const Edge& edge) const;
  Vertex target(const Edge& edge) const;

  /*! @brief Number of substituents of a vertex
   *
   * @complexity{@math{\Theta(1)}}
   */
  Vertex degree(Vertex a) const;

  //! Number of vertices in the graph
  Vertex N() const;
  //! Number of edges in the graph
  Vertex B() const;

  /*! @brief Checks whether all edges present in *this are present in @p other
   *
   * @complexity{@math{O(B)}}
   */
  bool identicalGraph(const PrivateGraph& other) const;

  /*! @brief Determine which vertices belong to which side of a bridge edge
   *
   * @complexity{@math{\Theta(N)}}
   * @note This function is not thread-safe.
   */
  std::pair<
    std::vector<AtomIndex>,
    std::vector<AtomIndex>
  > splitAlongBridge(Edge bridge) const;

//!@}

/*!
 * @name Cached properties access
 * Call complexity depends on whether the properties have been calculated
 * before. After generation, properties are valid as long as no non-const
 * method is called on this class that modifies state.
 *
 * None of these methods are thread-safe.
 * @{
 */
  void populateProperties() const;
  //! Access cycle information of the graph
  const Cycles& cycles() const;
  //! Access removal safety information of the graph
  const RemovalSafetyData& removalSafetyData() const;
  //! Access cycle information of the graph with eta bonds preserved
  const Cycles& etaPreservedCycles() const;
//!@}

//!@name Ranges
//!@{
  VertexRange vertices() const;
  EdgeRange edges() const;
  AdjacentVertexRange adjacents(Vertex a) const;
  IncidentEdgeRange edges(Vertex a) const;
//!@}

private:
//!@name Private types
//!@{
  struct Properties {
    boost::optional<RemovalSafetyData> removalSafetyDataOption;
    boost::optional<Cycles> cyclesOption;
    boost::optional<Cycles> etaPreservedCyclesOption;

    inline void invalidate() {
      removalSafetyDataOption = boost::none;
      cyclesOption = boost::none;
      etaPreservedCyclesOption = boost::none;
    }
  };

  RemovalSafetyData _generateRemovalSafetyData() const;
  Cycles _generateCycles() const;
  Cycles _generateEtaPreservedCycles() const;
//!@}

//!@name Information
//!@{
//!@}

//!@name Private state
//!@{
  //! A directly owned Boost Library Graph.
  BglType _graph;
  //! Property caching
  mutable Properties _properties;
//!@}
};

} // namespace molassembler
} // namespace Scine
#endif
