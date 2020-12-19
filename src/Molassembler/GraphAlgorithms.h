/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Algorithms to help in the interpretation of the molecular graph
 */

#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_ALGORITHMS_H

#include "Utils/Geometry/ElementTypes.h"
#include "Types.h"
#include "boost/functional/hash.hpp"

#include <unordered_map>
#include <vector>
#include <memory>
#include <limits>

namespace Scine {
namespace Molassembler {

// Forward-declarations
class Graph;

/*! @brief Calculates the graph distance from a single atom index to all others
 *
 * Performs a BFS through the entire graph starting at the supplied index and
 * records the distance of vertices it encounters.
 *
 * @complexity{@math{O(N)}}
 *
 * @throws std::out_of_range If source >= N()
 *
 * @returns A vector containing the distances of all vertices to the supplied
 *   index
 */
MASM_EXPORT std::vector<unsigned> distance(AtomIndex source, const Graph& graph);

struct MASM_EXPORT PredecessorMap {
  std::vector<AtomIndex> predecessors;

  /**
   * @brief Generate path from source to target vertex
   *
   * @param target Target vertex of shortest path
   *
   * @return Path starting at source and ending at target, including both
   *   source and target vertices
   */
  std::vector<AtomIndex> path(AtomIndex target) const;
};

/**
 * @brief Generates shortest paths to each vertex in a graph
 *
 * @param source Vertex to start at
 * @param graph Graph containing the vertex
 *
 * @throws std::out_of_range If source >= N()
 *
 * @return A flat predecessor map
 */
MASM_EXPORT PredecessorMap shortestPaths(AtomIndex source, const Graph& graph);

/**
 * @brief Base class for edit cost functors for graph edit distance methods
 *
 * This class defines the interface that a cost functor for the family of graph
 * edit distance methods must implement. Regarding the used terminology:
 * Alterations are insertion or deletion. Substitutions are changes to existing
 * parts of the graph.
 */
struct MASM_EXPORT EditCost {
  virtual ~EditCost() = default;

  // Alteration means insertion or deletion
  virtual unsigned vertexAlteration() const = 0;
  virtual unsigned edgeAlteration() const = 0;
  virtual unsigned elementSubstitution(Utils::ElementType from, Utils::ElementType to) const = 0;
  virtual unsigned bondSubstitution(BondType from, BondType to) const = 0;
};

/**
 * @brief Fuzzy-matching cost function that prefers element substitutions to
 *   alterations
 *
 * Graph edit distances generated with this cost function will prefer
 * substitutions to alterations, but not strongly.
 */
struct MASM_EXPORT FuzzyCost : public EditCost {
  unsigned vertexAlteration() const override { return 1; }
  unsigned edgeAlteration() const override { return 1; }
  unsigned elementSubstitution(Utils::ElementType from, Utils::ElementType to) const override {
    return (from != to) ? 1 : 0;
  }
  unsigned bondSubstitution(BondType from, BondType to) const override {
    return (from != to) ? 1 : 0;
  }
};

/**
 * @brief Element type conserving cost function, really allowing only edge
 *   alterations
 *
 * This cost function conserves element types and prevents vertex deletion.
 * Perfect for chemical reactions.
 */
struct MASM_EXPORT ElementsConservedCost : public FuzzyCost {
  unsigned vertexAlteration() const override { return 100; }
  unsigned elementSubstitution(Utils::ElementType from, Utils::ElementType to) const override {
    return (from != to) ? 100 : 0;
  }
};

/**
 * @brief Result type for single-pair graph edit distance calculation
 */
struct MASM_EXPORT MinimalGraphEdits {
  //! Flat vertex index map between the graphs
  using IndexMap = std::vector<AtomIndex>;

  //! Type for non-zero cost vertex edits in the result set
  struct VertexEdit {
    VertexEdit() = default;
    VertexEdit(AtomIndex i, AtomIndex j, unsigned c) : first(i), second(j), cost(c) {}

    //! Left graph index (possibly epsilon)
    AtomIndex first;
    //! Right graph index (possibly epsilon)
    AtomIndex second;
    //! Cost of the edit per the cost function
    unsigned cost;
  };

  //! Type for non-zero cost edge edits in the result set
  struct EdgeEdit {
    EdgeEdit() = default;
    EdgeEdit(BondIndex i, BondIndex j, unsigned c) : first(i), second(j), cost(c) {}

    //! Left graph bond index (possibly containing epsilon values)
    BondIndex first;
    //! Right graph bond index (possibly containing epsilon values)
    BondIndex second;
    //! Cost of the edit per the cost function
    unsigned cost;
  };

  //! Sentinel value for non-existent vertex
  static constexpr AtomIndex epsilon = std::numeric_limits<AtomIndex>::max();

  //! Total distance between graphs according to cost function
  unsigned distance;
  //! Vertex mapping, possibly including epsilon values
  IndexMap indexMap;
  //! List of non-zero-cost vertex edits
  std::vector<VertexEdit> vertexEdits;
  //! List of non-zero-cost edge edits
  std::vector<EdgeEdit> edgeEdits;
};

/**
 * @brief Graph edit distance calculation
 *
 * A maximum common connected subgraph is used to precondition the exact graph
 * edit distance algorithm. Otherwise, the exact graph edit distance quickly
 * becomes intractable past 10 vertices due to combinatorial space explosion
 * and rapid memory exhaustion.
 *
 * @param a First graph to calculate edit distance for
 * @param b Second graph to calculate edit distance for
 * @param cost Cost functor for graph edits
 *
 * @returns distance, index mapping and non-zero cost edit lists
 */
MASM_EXPORT MinimalGraphEdits minimalEdits(const Graph& a, const Graph& b, const EditCost& cost = FuzzyCost {});

//! Aggregate for multiple-graph graph edit distance
struct MASM_EXPORT MinimalReactionEdits {
  //! Input component and atom index therein
  using ComponentIndexPair = std::pair<unsigned, AtomIndex>;

  //! Type for non-zero cost vertex edits in the result set
  struct VertexEdit {
    VertexEdit() = default;
    VertexEdit(ComponentIndexPair i, ComponentIndexPair j, unsigned c)
      : first(std::move(i)), second(std::move(j)), cost(c) {}

    //! Left graph index and component
    ComponentIndexPair first;
    //! Right graph index and component
    ComponentIndexPair second;
    //! Cost of the edit per the cost function
    unsigned cost;
  };

  //! Type for non-zero cost edge edits in the result set
  struct EdgeEdit {
    //! Type for two component and index pairs
    using ComponentBondIndex = std::pair<ComponentIndexPair, ComponentIndexPair>;

    EdgeEdit() = default;
    EdgeEdit(ComponentBondIndex i, ComponentBondIndex j, unsigned c)
      : first(std::move(i)), second(std::move(j)), cost(c) {}

    //! Left bond
    ComponentBondIndex first;
    //! Right bond
    ComponentBondIndex second;
    //! Cost of the edit per the cost function
    unsigned cost;
  };

  //! Minimal edit cost between graphs
  unsigned distance;
  //! Vertex and component mapping
  std::unordered_map<ComponentIndexPair, ComponentIndexPair, boost::hash<ComponentIndexPair>> indexMap;
  //! List of non-zero-cost vertex edits
  std::vector<VertexEdit> vertexEdits;
  //! List of non-zero-cost edge edits
  std::vector<EdgeEdit> edgeEdits;
};

//! Input type for multiple graphs
using GraphList = std::vector<std::reference_wrapper<Graph>>;

/**
 * @brief Graph edit distance calculation for reactions
 *
 * A maximum common subgraph is used to precondition the exact graph edit
 * distance algorithm. Otherwise, the exact graph edit distance quickly becomes
 * intractable past 10 vertices due to combinatorial space explosion and rapid
 * memory exhaustion.
 *
 * @param lhs List of graphs of the left side of the reaction
 * @param rhs List of graphs of the right side of the reaction
 *
 * @throws std::logic_error If the number of atoms in both sides is unequal or
 * the element composition of both sides is different.
 *
 * @returns distance, index mapping and non-zero cost edit lists
 */
MASM_EXPORT MinimalReactionEdits reactionEdits(const GraphList& lhsGraphs, const GraphList& rhsGraphs);

} // namespace Molassembler
} // namespace Scine

#endif
