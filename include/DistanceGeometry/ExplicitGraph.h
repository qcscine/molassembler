#ifndef INCLUDE_DG_LIMITS_GRAPH_H
#define INCLUDE_DG_LIMITS_GRAPH_H

#include "Molecule.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"

/*! @file
 *
 * Declaration of a graph class that aids in the determination of triangle
 * inequality limit distance bounds and the generation of a distance matrix
 * from a list of atom-index pairwise distance bounds.
 */

namespace MoleculeManip {

namespace DistanceGeometry {

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
  using ExplicitGraphType = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::directedS,
    boost::no_property,
    EdgeWeightProperty
  >;

private:
  ExplicitGraphType _graph;
  bool _generatedDistances;
  const Molecule& _molecule;

  void _updateOrAddEdge(
    const ExplicitGraphType::vertex_descriptor& a,
    const ExplicitGraphType::vertex_descriptor& b,
    const double& edgeWeight
  );

  void _updateGraphWithFixedDistance(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const double& fixedDistance
  );

  static inline ExplicitGraphType::vertex_descriptor left(const AtomIndexType& a) {
    return 2 * a;
  }

  static inline ExplicitGraphType::vertex_descriptor right(const AtomIndexType& a) {
    return 2 * a + 1;
  }

public:
  using BoundList = std::vector<
    std::tuple<AtomIndexType, AtomIndexType, ValueBounds>
  >;
  ExplicitGraph(const Molecule& molecule, const BoundList& bounds);

  //! Adds edges to the underlying graph to represent the bound between the atoms
  void addBound(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const ValueBounds& bound
  );

  //! Adds edges to the underlying graph representing implicit lower bounds
  void addImplicitEdges();

  void setDistance(
    const ExplicitGraphType::vertex_descriptor& a,
    const ExplicitGraphType::vertex_descriptor& b,
    double distance
  );

  double lowerBound(
    const AtomIndexType& i,
    const AtomIndexType& j
  ) const;

  double upperBound(
    const AtomIndexType& i,
    const AtomIndexType& j
  ) const;


  const ExplicitGraphType& getGraph() const;

  DistanceBoundsMatrix makeDistanceBounds() const;

  /*!
   * Generates a distances matrix conforming to the triangle inequality bounds
   * while modifying state information. Can only be called once!
   */
  Eigen::MatrixXd makeDistanceMatrix();
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
