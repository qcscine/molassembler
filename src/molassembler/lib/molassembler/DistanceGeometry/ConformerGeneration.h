// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_CONFORMER_GENERATION_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_CONFORMER_GENERATION_H

#include "boost_outcome/outcome.hpp"

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
  const Eigen::VectorXd& vectorizedPositions
);

AngstromWrapper convertToAngstromWrapper(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
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

/*!
 * A logging, not throwing, otherwise identical implementation of
 * runDistanceGeometry that returns detailed intermediate data from
 * refinements, while run returns only the final conformers.
 */
std::list<RefinementData> debug(
  const Molecule& molecule,
  unsigned numConformers,
  Partiality metrizationOption = Partiality::FourAtom,
  bool useYInversionTrick = true
);

/*!
 * The main implementation of Distance Geometry. Generates an ensemble of 3D
 * structures of a given Molecule.
 *
 * @param metrizationOption After choosing an element of distance matrix between its
 * triangle inequality bounds, it is optional whether to ensure that all other
 * bounds afterwards also conform to the triangle inequality. Since the slack
 * removed from the distance bounds per chosen distance and thus the accuracy
 * gained decrease exponentially, you may choose to perform re-smoothing only
 * for a limited set of atoms.
 *
 * @param useYInversionTrick After embedding coordinates for the first time,
 * whether chiral constraints are correct by sign is normally distributed
 * around 0.5. If fewer than half of all chiral constraints are correct, an
 * inversion of a coordinate will lead to a structure that has exactly 1 - x
 * chiral constraints correct.
 */
outcome::result<
  std::vector<AngstromWrapper>
> run(
  const Molecule& molecule,
  unsigned numConformers,
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
