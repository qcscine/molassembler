// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_CONFORMER_GENERATION_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_CONFORMER_GENERATION_H

#include "boost_outcome/outcome.hpp"

#include "molassembler/Conformers.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/RefinementDebugData.h"
#include "molassembler/Log.h"

/*! @file
 *
 * @brief Distance Geometry conformation generating procedures
 *
 * Declares the central conformation (and -ensemble) generating functions that
 * start the DG procedure.
 */

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

namespace detail {

AngstromWrapper convertToAngstromWrapper(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
);

Eigen::MatrixXd fitAndSetFixedPositions(
  const dlib::matrix<double, 0, 1>& vectorizedPositions,
  const Configuration& configuration
);

/*!
 * @brief Assigns any unassigned stereopermutators in a molecule at random
 *
 * If any stereopermutators in a molecule are unassigned, we must progressively
 * assign them at random (consistent with relative occurrences) through
 * Molecule's interface (so that ranking change effects w.r.t. the number of
 * stereopermutations in permutators are handled gracefully) before attempting
 * to model the Molecule (which requires that all stereopermutators are
 * assigned).
 */
Molecule narrow(Molecule molecule);

} // namespace detail

//! Intermediate conformational data about a Molecule given by a spatial model
struct MoleculeDGInformation {
  SpatialModel::BoundsList bounds;
  std::vector<ChiralityConstraint> chiralityConstraints;
};

//! Collects intermediate conformational data about a Molecule using a spatial model
MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const Configuration& configuration,
  double looseningFactor = 1.0
);

//! Debug function, also collects graphviz of the conformational model
MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const Configuration& configuration,
  double looseningFactor,
  std::string& spatialModelGraphvizString
);

/*!
 * A logging, not throwing, mostly identical implementation of
 * runDistanceGeometry that returns detailed intermediate data from
 * refinements, while run returns only the final conformers.
 */
std::list<RefinementData> debug(
  const Molecule& molecule,
  unsigned numConformers,
  const Configuration& configuration
);

/*!
 * The main implementation of Distance Geometry. Generates an ensemble of 3D
 * structures of a given Molecule.
 *
 * @param configuration A configuration object controlling the Distance
 *   Geometry procedure
 */
outcome::result<
  std::vector<AngstromWrapper>
> run(
  const Molecule& molecule,
  unsigned numConformers,
  const Configuration& configuration
);

} // namespace DistanceGeometry

} // namespace molassembler

#endif
