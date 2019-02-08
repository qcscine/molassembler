/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Boost Graph Library wrapper to help in concealing underlying type
 */

#ifndef INCLUDE_MOLASSEMBLER_INNER_GRAPH_H
#define INCLUDE_MOLASSEMBLER_INNER_GRAPH_H

#include "boost/optional/optional_fwd.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "Utils/ElementTypes.h"

#include "molassembler/Cycles.h"
#include "molassembler/Types.h"

#include <limits>

namespace Scine {

namespace molassembler {

/**
 * @brief Library internal graph class wrapping BGL types
 */
class InnerGraph {
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
  using BGLType = boost::adjacency_list<
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

  using Vertex = BGLType::vertex_descriptor;
  using Edge = BGLType::edge_descriptor;

  template<typename Iter>
  using Range = std::pair<Iter, Iter>;

  using VertexRange = Range<BGLType::vertex_iterator>;
  using EdgeRange = Range<BGLType::edge_iterator>;
  using AdjacentVertexRange = Range<BGLType::adjacency_iterator>;
  using IncidentEdgeRange = Range<BGLType::out_edge_iterator>;

  /*!
   * @brief Data class to return removal safety information on the graph
   *
   * @note These are grouped because it is fairly easy to work out bridges from
   * articulation vertices and hence their joint calculation makes a lot of
   * sense. It might also be easier to cache this way without the dependence.
   */
  struct RemovalSafetyData {
    //! Articulation vertices cannot be removed without disconnecting the graph
    std::unordered_set<InnerGraph::Vertex> articulationVertices;
    //! Bridges are edges that cannot be removed without disconnecting the graph
    std::set<InnerGraph::Edge> bridges;
  };
//!@}

//!@name Constructors
//!@{
  InnerGraph();
  InnerGraph(Vertex N);
//!@}

//!@name Static members
//!@{
  static constexpr Vertex removalPlaceholder = std::numeric_limits<Vertex>::max();
//!@}

//!@name Modification
//!@{
  Edge addEdge(Vertex a, Vertex b, BondType bondType);

  Vertex addVertex(Utils::ElementType elementType);

  void applyPermutation(const std::vector<Vertex>& permutation);

  BondType& bondType(const Edge& edge);

  BGLType& bgl();

  //! Removes all bonds involving a vertex
  void clearVertex(Vertex a);

  Utils::ElementType& elementType(Vertex a);

  //! Removes an edge from the graph.
  void removeEdge(const Edge& e);

  /*! Removes a vertex from the graph.
   *
   * @warning This invalidates ALL vertex and edge descriptors!
   */
  void removeVertex(Vertex a);
//!@}

//!@name Information
//!@{
  Utils::ElementType elementType(Vertex a) const;
  BondType bondType(const Edge& edge) const;
  const BGLType& bgl() const;

  /*!
   * @brief Determine whether a vertex can be safely removed
   * A Vertex can be safely removed if it is not an articulation vertex
   * @note This function is not thread-safe.
   */
  bool canRemove(Vertex a) const;
  /*!
   * @brief Determine whether an edge can be safely removed
   * An edge can be safely removed if it is not a bridge edge
   * @note This function is not thread-safe.
   */
  bool canRemove(const Edge& edge) const;

  unsigned connectedComponents() const;
  unsigned connectedComponents(std::vector<unsigned>& componentMap) const;

  //! Precondition: Edge exists!
  Edge edge(Vertex a, Vertex b) const;
  boost::optional<Edge> edgeOption(Vertex a, Vertex b) const;
  Vertex source(const Edge& edge) const;
  Vertex target(const Edge& edge) const;
  Vertex degree(Vertex a) const;
  Vertex N() const;
  Vertex B() const;

  //! Checks whether all edges present in *this are present in @p other
  bool identicalGraph(const InnerGraph& other) const;

  /*!
   * @brief Determine which vertices belong to which side of a bridge edge
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
  //! Access cycle information of the graph
  const Cycles& cycles() const;
  //! Access removal safety information of the graph
  const RemovalSafetyData& removalSafetyData() const;
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

    inline void invalidate() {
      removalSafetyDataOption = boost::none;
      cyclesOption = boost::none;
    }
  };
//!@}

//!@name Information
//!@{
//!@}

//!@name Private state
//!@{
  //! A directly owned Boost Library Graph.
  BGLType _graph;
  //! Property caching
  mutable Properties _properties;
//!@}
};

} // namespace molassembler

} // namespace Scine
#endif
