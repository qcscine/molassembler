/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Fast shortest paths graph for triangle inequalities
 *
 * Declaration of a graph class that aids in the determination of triangle
 * inequality limit distance bounds and the generation of a distance matrix
 * from a list of atom-index pairwise distance bounds.
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_EXPLICIT_GRAPH_H
#define INCLUDE_MOLASSEMBLER_DG_EXPLICIT_GRAPH_H

// #define MOLASSEMBLER_EXPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM

#include "boost/graph/adjacency_list.hpp"
#include "boost_outcome/outcome.hpp"
#include "Eigen/Core"
#include "Utils/ElementInfo.h"
#include "temple/Preprocessor.h"

#include "molassembler/DistanceGeometry/ValueBounds.h"

namespace Scine {

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

// Forward-declare Molecule
class Molecule;

namespace DistanceGeometry {

// Forward-declarations
enum class Partiality;
class DistanceBoundsMatrix;


/*!
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
class ExplicitGraph {
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

  using BoundsList = std::map<
    std::array<VertexDescriptor, 2>,
    ValueBounds
  >;
//!@}

//!@name Special member functions
//!@{
  ExplicitGraph(
    const Molecule& molecule,
    const DistanceBoundsMatrix& bounds
  );

  ExplicitGraph(
    const Molecule& molecule,
    const BoundsList& bounds
  );
//!@}

//!@name Static member functions
//!@{
  static void _explainContradictionPaths(
    VertexDescriptor a,
    VertexDescriptor b,
    const std::vector<VertexDescriptor>& predecessors,
    const std::vector<double>& distances
  );

  PURITY_STRONG static inline bool isLeft(const VertexDescriptor i) {
    return i % 2 == 0;
  }

  PURITY_STRONG static inline VertexDescriptor left(const VertexDescriptor a) {
    return 2 * a;
  }

  PURITY_STRONG static inline VertexDescriptor right(const VertexDescriptor a) {
    return 2 * a + 1;
  }
//!@}

//!@name Modifiers
//!@{
  void addBound(
    VertexDescriptor a,
    VertexDescriptor b,
    const ValueBounds& bound
  );

  //! Adds edges to the underlying graph representing implicit lower bounds
  void addImplicitEdges();

  void setDistance(
    VertexDescriptor a,
    VertexDescriptor b,
    double distance
  );

  double lowerBound(
    VertexDescriptor a,
    VertexDescriptor b
  ) const;

  double upperBound(
    VertexDescriptor a,
    VertexDescriptor b
  ) const;

  /*!
   * Generates a distances matrix conforming to the triangle inequality bounds
   * while modifying state information. Can only be called once!
   */
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix() noexcept;

  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(Partiality partiality) noexcept;
//!@}

//!@name Information
//!@{
  //! Returns the length of the maximal implicit lower bound outgoing from a left vertex
  double maximalImplicitLowerBound(VertexDescriptor i) const;

  const GraphType& graph() const;

  outcome::result<Eigen::MatrixXd> makeDistanceBounds() const noexcept;
//!@}

private:
  GraphType _graph;
  const Molecule& _molecule;
  //! Stores the two heaviest element types
  std::array<Scine::Utils::ElementType, 2> _heaviestAtoms;

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

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine

#endif
