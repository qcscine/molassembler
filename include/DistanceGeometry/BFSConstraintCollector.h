#ifndef INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H
#define INCLUDE_DG_BFS_CONSTRAINT_COLLECTOR_H

#include "DistanceGeometry/DistanceBoundsMatrix.h"
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

  enum class DistanceMethod {
    Uniform,
    UFFLike
  };

private:
  /* Private Members */
  // Pre-set constants
  const double oneTwoVariance = 0.1;

  // Closures
  const AdjacencyList& _adjacencies;
  const DistanceMethod& _distanceMethod;

  // Output (via function operator side effect)
  DistanceBoundsMatrix& _distanceBounds;

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
  ) {
    return _stereocenterMap[j] -> angle(i, j, k);
  }

public:
  /* Ctor */
  /*! 
   * Assumes the Distance Bounds is entirely uninitialized (lower diagonal is 0,
   * upper diagonal is 100).
   */
  BFSConstraintCollector(
    const AdjacencyList& adjacencies,
    const StereocenterList& stereocenterList,
    DistanceBoundsMatrix& distanceBounds,
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

  /*! After all trees have been visited, replace all lower bounds that are 0
   * to the sum of the respective atoms' vdw radii and then smooth the matrix.
   * Also checks that no bounds inconsistencies exists afterwards in debug
   * builds.
   */
  void finalizeBoundsMatrix();

/* Information */
  // Collect chirality prototypes from gathered information
  std::vector<
    Stereocenters::ChiralityConstraintPrototype
  > getChiralityPrototypes() const;
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
