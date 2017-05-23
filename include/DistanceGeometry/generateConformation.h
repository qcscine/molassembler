#ifndef INCLUDE_DG_GENERATE_CONFORMATION_H
#define INCLUDE_DG_GENERATE_CONFORMATION_H

#include "Molecule.h"
#include "Delib/PositionCollection.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/BFSConstraintCollector.h"
#include "DistanceGeometry/RefinementDebugData.h"
#include "Log.h"

#include <vector>
#include <Eigen/Core>

/* TODO
 */

namespace MoleculeManip {

namespace DistanceGeometry {

namespace detail {

Delib::PositionCollection convertToPositionCollection(
  const Eigen::VectorXd& vectorizedPositions
);

Delib::PositionCollection convertToPositionCollection(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
);

ChiralityConstraint propagate(
  const DistanceBoundsMatrix& bounds,
  const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
);

DGDebugData debugDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const bool& useYInversionTrick = true,
  const BFSConstraintCollector::DistanceMethod& distanceMethod = BFSConstraintCollector::DistanceMethod::UFFLike
);

std::list<Delib::PositionCollection> runDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const bool& useYInversionTrick = true,
  const BFSConstraintCollector::DistanceMethod& distanceMethod = BFSConstraintCollector::DistanceMethod::UFFLike
);

} // namespace detail

struct MoleculeDGInformation {
  DistanceBoundsMatrix distanceBounds;
  std::vector<
    Stereocenters::Stereocenter::ChiralityConstraintPrototype
  > chiralityConstraintPrototypes;

  explicit MoleculeDGInformation(const unsigned& N);
};

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const BFSConstraintCollector::DistanceMethod& distanceMethod = BFSConstraintCollector::DistanceMethod::UFFLike
);

// Public functions
std::list<Delib::PositionCollection> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization = MetrizationOption::off
);

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization = MetrizationOption::off
);

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
