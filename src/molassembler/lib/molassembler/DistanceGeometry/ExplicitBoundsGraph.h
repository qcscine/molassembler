/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Fast shortest paths graph for triangle inequalities
 *
 * Declaration of a graph class that aids in the determination of triangle
 * inequality limit distance bounds and the generation of a distance matrix
 * from a list of atom-index pairwise distance bounds.
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_EXPLICIT_BOUNDS_GRAPH_H
#define INCLUDE_MOLASSEMBLER_DG_EXPLICIT_BOUNDS_GRAPH_H

// #define MOLASSEMBLER_EXPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM

#include "boost/graph/adjacency_list.hpp"
#include "boost_outcome/outcome.hpp"
#include "Eigen/Core"
#include "Utils/Geometry/ElementInfo.h"

#include "molassembler/Export.h"
#include "molassembler/DistanceGeometry/ValueBounds.h"

namespace Scine {

namespace molassembler {

namespace random {
class Engine;
} // namespace random

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

// Forward-declare PrivateGraph
class PrivateGraph;

namespace distance_geometry {

// Forward-declarations
enum class Partiality;
class DistanceBoundsMatrix;


/*! @brief BGL wrapper to help with distance bounds smoothing
 *
 * A class that helps with lower complexity determination of distance bounds
 * compatible with the triangle inequality limits via reinterpretation as a
 * graph shortest paths calculation.
 *
 * The main idea is to transfer all ValueBounds collected by a
 * SpatialModel to this class. The compatible triangle inequality bounds
 * to this set of ValueBounds should then be extracted as a DistanceBoundsMatrix
 * and kept for refinement. For every new conformation, the underlying graph
 * containing the information from the ValueBounds needs to be copied and
 * generateDistanceMatrix called upon it. This procedure modifies the
 * underlying graph and hence cannot be called repeatedly.
 *
 * The underlying data structure is a fully explicit BGL graph containing all
 * edges and edge weights.
 */
class ExplicitBoundsGraph {
public:
//!@name Member types
//!@{
  using EdgeWeightProperty = boost::property<boost::edge_weight_t, double>;
  using GraphType = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::directedS,
    boost::no_property,
    EdgeWeightProperty
  >;
  using VertexDescriptor = GraphType::vertex_descriptor;
  using EdgeDescriptor = GraphType::edge_descriptor;

  using BoundsMatrix = Eigen::MatrixXd;
//!@}

//!@name Special member functions
//!@{
  /*! @brief Construct from distance bounds matrix
   *
   * @complexity{@math{\Theta(N^2)}}
   */
  ExplicitBoundsGraph(
    const PrivateGraph& inner,
    const DistanceBoundsMatrix& bounds
  );

  /*! @brief Construct from bounds matrix
   *
   * @complexity{@math{\Theta(N^2)}}
   */
  ExplicitBoundsGraph(
    const PrivateGraph& inner,
    const BoundsMatrix& bounds
  );
//!@}

//!@name Static member functions
//!@{
  /*! @brief Explain contradictions between distance bounds on a path
   *
   * @complexity{linear in the length of the path between @p a and @p b}
   */
  static void explainContradictionPaths(
    VertexDescriptor a,
    VertexDescriptor b,
    const std::vector<VertexDescriptor>& predecessors,
    const std::vector<double>& distances
  );

  /*! @brief Check whether a vertex is part of the left subgraph
   *
   * @complexity{@math{\Theta(1)}}
   */
  [[gnu::const]] static inline bool isLeft(const VertexDescriptor i) noexcept {
    return i % 2 == 0;
  }

  [[gnu::const]] static inline bool sameSide(const VertexDescriptor i, const VertexDescriptor j) noexcept {
    return i % 2 == j % 2;
  }

  /*! @brief Get the left subgraph vertex descriptor corresponding to an outer index
   *
   * @complexity{@math{\Theta(1)}}
   */
  [[gnu::const]] static inline VertexDescriptor left(const VertexDescriptor a) noexcept {
    return 2 * a;
  }

  /*! @brief Get the right subgraph vertex descriptor corresponding to an outer index
   *
   * @complexity{@math{\Theta(1)}}
   */
  [[gnu::const]] static inline VertexDescriptor right(const VertexDescriptor a) noexcept {
    return 2 * a + 1;
  }
//!@}

//!@name Modifiers
//!@{
  /*! @brief Adds a bound between outer vertex indices to the graph
   *
   * This is represented by six edges:
   * - Bidirectional within the left subgraph between a and b weighted with the
   *   upper bound
   * - Bidirectional within the right subgraph between a and b weighted with
   *   the upper bound
   * - Edges from the left vertices of a and b to the right opposite one with
   *   the negative lower bound
   *
   * @complexity{@math{\Theta(1)}}
   */
  void addBound(
    VertexDescriptor a,
    VertexDescriptor b,
    const ValueBounds& bound
  );

  /*! @brief Fetches the lower bound between outer vertex indices from the graph
   *
   * @complexity{@math{\Theta(1)}}
   */
  double lowerBound(VertexDescriptor a, VertexDescriptor b) const;

  /*! @brief Fetches the upper bound between outer vertex indices from the graph
   *
   * @complexity{@math{\Theta(1)}}
   */
  double upperBound(VertexDescriptor a, VertexDescriptor b) const;

  /*! @brief Generate a distance matrix
   *
   * Generates a distances matrix conforming to the triangle inequality bounds
   * while modifying state information. Can only be called once!
   *
   * @complexity{@math{O(V^2 \cdot E)}}
   */
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(random::Engine& engine) noexcept;

  //!@overload
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(random::Engine& engine, Partiality partiality) noexcept;
//!@}

//!@name Information
//!@{
  /*! @brief Returns the length of the maximal implicit lower bound outgoing from a left vertex
   *
   * @complexity{@math{\Theta(1)}}
   */
  double maximalImplicitLowerBound(VertexDescriptor i) const;

  /*! @brief Nonmodifiable access to underlying graph
   *
   * @complexity{@math{\Theta(1)}}
   */
  const GraphType& graph() const;

  /*! @brief Make smooth distance bounds
   *
   * @complexity{@math{\Theta(V \cdot E)}}
   */
  outcome::result<Eigen::MatrixXd> makeDistanceBounds() const noexcept;
//!@}

private:
  GraphType _graph;
  const PrivateGraph& _inner;
  //! Stores the two heaviest element types
  std::array<Utils::ElementType, 2> _heaviestAtoms;

  void _updateOrAddEdge(
    VertexDescriptor i,
    VertexDescriptor j,
    double edgeWeight
  );

  void _updateGraphWithFixedDistance(
    VertexDescriptor a,
    VertexDescriptor b,
    double fixedDistance
  );
};

} // namespace distance_geometry

} // namespace molassembler

} // namespace Scine

#endif
