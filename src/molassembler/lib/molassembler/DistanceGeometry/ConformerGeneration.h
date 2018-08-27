#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_CONFORMER_GENERATION_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_CONFORMER_GENERATION_H

#include "boost_outcome/outcome.hpp"

#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/RefinementDebugData.h"
#include "molassembler/Log.h"

/*! @file
 *
 * Declares the central conformation (and -ensemble) generating functions that
 * start the DG procedure.
 */

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

namespace detail {

AngstromWrapper convertToAngstromWrapper(
  const Eigen::VectorXd& vectorizedPositions
);

AngstromWrapper convertToAngstromWrapper(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
);

//! Assigns any unassigned stereocenters in a molecule at random
Molecule narrow(Molecule moleculeCopy);

} // namespace detail

/*!
 * A logging, not throwing otherwise identical implementation of
 * runDistanceGeometry, that returns detailed intermediate data from a
 * refinement, while runDistanceGeometry returns only the final result.
 *
 * @note Contained PositionCollections are in Angstrom length units
 */
std::list<RefinementData> debug(
  const Molecule& molecule,
  unsigned numStructures,
  Partiality metrizationOption = Partiality::FourAtom,
  bool useYInversionTrick = true
);

/*!
 * The main implementation of Distance Geometry. Generates an ensemble of 3D
 * structures of a given Molecule.
 *
 * Metrization options: After choosing an element of distance matrix between its
 * triangle inequality bounds, it is optional whether to ensure that all other
 * bounds afterwards also conform to the triangle inequality. Since the slack
 * removed from the distance bounds per chosen distance and thus the accuracy
 * gained decrease exponentially, you may choose to perform re-smoothing only
 * for a limited set of atoms.
 *
 * Use Y-inversion trick: After embedding coordinates for the first time,
 * whether chiral constraints are correct by sign is normally distributed
 * around 0.5. If fewer than half of all chiral constraints are correct, an
 * inversion of a coordinate will lead to a structure that has exactly 1 - x
 * chiral constraints correct.
 *
 * Distance method: For debug purposes, using uniform distances between atoms
 * may be desirable for particularly hypothetical structures.
 *
 * @note Returns PositionCollections in Angstrom length units
 */
outcome::result<
  std::vector<AngstromWrapper>
> run(
  const Molecule& molecule,
  unsigned numStructures,
  Partiality metrizationOption = Partiality::FourAtom,
  bool useYInversionTrick = true
);

//! Intermediate conformational data about a Molecule given by a spatial model
struct MoleculeDGInformation {
  SpatialModel::BoundsList bounds;
  std::vector<ChiralityConstraint> chiralityConstraints;
};

//! Collects intermediate conformational data about a Molecule using a spatial model
MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  double looseningFactor = 1.0
);

//! Debug function, also collects graphviz of the conformational model
MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  double looseningFactor,
  std::string& spatialModelGraphvizString
);

} // namespace DistanceGeometry

} // namespace molassembler

#endif
