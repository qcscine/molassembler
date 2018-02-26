#ifndef INCLUDE_DG_LIMITS_GRAPH_H
#define INCLUDE_DG_LIMITS_GRAPH_H

#include "DistanceGeometry/ValueBounds.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/outcome.hpp"
#include "Eigen/Core"
#include "Delib/ElementInfo.h"

/*! @file
 *
 * Declaration of a graph class that aids in the determination of triangle
 * inequality limit distance bounds and the generation of a distance matrix
 * from a list of atom-index pairwise distance bounds.
 */

namespace MoleculeManip {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

// Forward-declare Molecule
class Molecule;

namespace DistanceGeometry {

// Forward-declare Partiality
enum class Partiality;

/*!
 * A class that helps with lower complexity determination of distance bounds
 * compatible with the triangle inequality limits via reinterpretation as a
 * graph shortest paths calculation.
 *
 * The main idea is to transfer all ValueBounds collected by a
 * MoleculeSpatialModel to this class. The compatible triangle inequality bounds
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

private:
  GraphType _graph;
  const Molecule& _molecule;
  //! Stores the two heaviest element types
  std::array<Delib::ElementType, 2> _heaviestAtoms;

  void _updateOrAddEdge(
    const VertexDescriptor& a,
    const VertexDescriptor& b,
    const double& edgeWeight
  );

  void _updateGraphWithFixedDistance(
    const VertexDescriptor& a,
    const VertexDescriptor& b,
    const double& fixedDistance
  );

  static inline VertexDescriptor left(const VertexDescriptor& a) {
    return 2 * a;
  }

  static inline VertexDescriptor right(const VertexDescriptor& a) {
    return 2 * a + 1;
  }

public:
  using BoundList = std::vector<
    std::tuple<VertexDescriptor, VertexDescriptor, ValueBounds>
  >;
  ExplicitGraph(const Molecule& molecule, const BoundList& bounds);

  static inline bool isLeft(const VertexDescriptor& i) {
    return i % 2 == 0;
  }

  //! Adds edges to the underlying graph to represent the bound between the atoms
  void addBound(
    const VertexDescriptor& a,
    const VertexDescriptor& b,
    const ValueBounds& bound
  );

  //! Adds edges to the underlying graph representing implicit lower bounds
  void addImplicitEdges();

  void setDistance(
    const VertexDescriptor& a,
    const VertexDescriptor& b,
    double distance
  );

  double lowerBound(
    const VertexDescriptor& i,
    const VertexDescriptor& j
  ) const;

  double upperBound(
    const VertexDescriptor& i,
    const VertexDescriptor& j
  ) const;

  //! Returns the length of the maximal implicit lower bound outgoing from a left vertex
  double maximalImplicitLowerBound(const VertexDescriptor& i) const;

  const GraphType& getGraph() const;

  outcome::result<Eigen::MatrixXd> makeDistanceBounds() const noexcept;

  /*!
   * Generates a distances matrix conforming to the triangle inequality bounds
   * while modifying state information. Can only be called once!
   */
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix() noexcept;

  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(Partiality partiality) noexcept;
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
