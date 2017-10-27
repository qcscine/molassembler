#ifndef INCLUDE_DG_GENERATE_CONFORMATION_H
#define INCLUDE_DG_GENERATE_CONFORMATION_H

#include "Delib/PositionCollection.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "DistanceGeometry/RefinementDebugData.h"
#include "Log.h"
#include "Molecule.h"

#include <vector>
#include <Eigen/Core>

/*! @file
 *
 * Declares the central conformation (and -ensemble) generating functions that
 * start the DG procedure.
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
  const Stereocenters::ChiralityConstraintPrototype& prototype
);

DGDebugData debugDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const bool& useYInversionTrick = true,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod = MoleculeSpatialModel::DistanceMethod::UFFLike
);

std::list<Delib::PositionCollection> runDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const bool& useYInversionTrick = true,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod = MoleculeSpatialModel::DistanceMethod::UFFLike
);

} // namespace detail

struct MoleculeDGInformation {
  DistanceBoundsMatrix distanceBounds;
  std::vector<
    Stereocenters::ChiralityConstraintPrototype
  > chiralityConstraintPrototypes;

  explicit MoleculeDGInformation(const unsigned& N);
};

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod = MoleculeSpatialModel::DistanceMethod::UFFLike
);

// "Public" functions
std::list<Delib::PositionCollection> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization = MetrizationOption::full
);

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization = MetrizationOption::full
);

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
