/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
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
#include "Eigen/Core"
#include "Utils/Geometry/ElementInfo.h"

#include "Molassembler/DistanceGeometry/ValueBounds.h"
#include "Molassembler/Conformers.h"

namespace Scine {
namespace Molassembler {

namespace Random {
class Engine;
} // namespace Random

// Forward-declare PrivateGraph
class PrivateGraph;

namespace DistanceGeometry {

// Forward-declarations
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
  Result<Eigen::MatrixXd> makeDistanceMatrix(Random::Engine& engine) noexcept;

  //!@overload
  Result<Eigen::MatrixXd> makeDistanceMatrix(Random::Engine& engine, Partiality partiality) noexcept;
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
  Result<Eigen::MatrixXd> makeDistanceBounds() const noexcept;
//!@}

private:
  GraphType graph_;
  const PrivateGraph& inner_;
  //! Stores the two heaviest element types
  std::array<Utils::ElementType, 2> heaviestAtoms_;

  void updateOrAddEdge_(
    VertexDescriptor i,
    VertexDescriptor j,
    double edgeWeight
  );

  void updateGraphWithFixedDistance_(
    VertexDescriptor a,
    VertexDescriptor b,
    double fixedDistance
  );
};

} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine

#endif
