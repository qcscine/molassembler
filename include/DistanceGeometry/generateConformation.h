#ifndef INCLUDE_DG_GENERATE_CONFORMATION_H
#define INCLUDE_DG_GENERATE_CONFORMATION_H

#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "DistanceGeometry/RefinementDebugData.h"
#include "Log.h"
#include "boost/outcome.hpp"

/*! @file
 *
 * Declares the central conformation (and -ensemble) generating functions that
 * start the DG procedure.
 */

namespace MoleculeManip {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

namespace detail {

Delib::PositionCollection convertToPositionCollection(
  const Eigen::VectorXd& vectorizedPositions
);

Delib::PositionCollection convertToPositionCollection(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
);

outcome::result<ChiralityConstraint> propagate(
  const DistanceBoundsMatrix& bounds,
  const Stereocenters::ChiralityConstraintPrototype& prototype
);

std::list<RefinementData> debugDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const Partiality& metrizationOption = Partiality::FourAtom,
  const bool& useYInversionTrick = true,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod = MoleculeSpatialModel::DistanceMethod::UFFLike
);

outcome::result<
  std::list<Delib::PositionCollection>
> runDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const Partiality& metrizationOption = Partiality::FourAtom,
  const bool& useYInversionTrick = true,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod = MoleculeSpatialModel::DistanceMethod::UFFLike
);

} // namespace detail


struct MoleculeDGInformation {
  MoleculeSpatialModel::BoundList boundList;
  std::vector<Stereocenters::ChiralityConstraintPrototype> chiralityConstraintPrototypes;
};

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod = MoleculeSpatialModel::DistanceMethod::UFFLike
);

// "Public" functions
outcome::result<
  std::list<Delib::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures
);

outcome::result<Delib::PositionCollection> generateConformation(const Molecule& molecule);

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
