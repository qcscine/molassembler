#ifndef INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H
#define INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H

#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "AdjacencyList.h"
#include "Tree.h"
#include "SequenceSet.h"

/* TODO
 * - tests
 */

namespace MoleculeManip {

namespace DistanceGeometry {

class BFSConstraintCollector {
public:
/* Typedefs */
  using NodeType = Tree::Node<AtomIndexType>;
  using SmallestCycleMapType = std::map<AtomIndexType, unsigned>;

  enum class DistanceMethod {
    Uniform,
    UFFLike
  };

private:
/* Private Members */
  // Closures
  const AdjacencyList& _adjacencies;
  const DistanceMethod& _distanceMethod;

  // Output (via function operator side effect)
  MoleculeSpatialModel& _spatialModel;

  // Constant members
  const CycleData _cycleData;
  const SmallestCycleMapType _smallestCycleMap;

  // State
  std::map<
    AtomIndexType,
    std::shared_ptr<Stereocenters::Stereocenter>
  > _stereocenterMap;

  SequenceSet<3> _angleSequences;
  SequenceSet<4> _dihedralSequences;

  inline double _getAngle(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) const {
    return _stereocenterMap.at(j) -> angle(i, j, k);
  }

  std::map<AtomIndexType, unsigned> _makeSmallestCycleMap() const;

  std::vector<AtomIndexType> _makeRingIndexSequence(
    const std::set<GraphType::edge_descriptor>& edgeSet
  ) const;

  unsigned _countPlanarityEnforcingBonds(
    const std::set<GraphType::edge_descriptor>& edgeSet
  ) const;

public:
  // Static constants
  /*! Relative bond distance variance, 0.0x meaning x% variance. Must fulfill
   * 0 < x << 1
   */
  static constexpr double bondRelativeVariance = 0.03; 
  //! Absolute angle variance in radians. Must fulfill 0 < x << M_PI
  static constexpr double angleAbsoluteVariance = M_PI / 36;

  /* Ctor */
  /*! 
   * Assumes the Distance Bounds is entirely uninitialized (lower diagonal is 0,
   * upper diagonal is 100).
   */
  BFSConstraintCollector(
    const AdjacencyList& adjacencies,
    const StereocenterList& stereocenterList,
    MoleculeSpatialModel& spatialModel,
    const DistanceMethod& distanceMethod = DistanceMethod::UFFLike
  );

/* Modifiers */
  //! Impure Function operator for Tree BFSVisit
  bool operator() (
    const std::shared_ptr<NodeType>& nodePtr,
    const unsigned& depth
  );

  void set12Bounds();
  void set13Bounds(const std::vector<AtomIndexType>& chain);
  void set14Bounds(const std::vector<AtomIndexType>& chain);

/* Information */
  // Collect chirality prototypes from gathered information
  std::vector<
    Stereocenters::ChiralityConstraintPrototype
  > getChiralityPrototypes() const;
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
