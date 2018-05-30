#ifndef INCLUDE_DG_GENERATE_CONFORMATION_H
#define INCLUDE_DG_GENERATE_CONFORMATION_H

#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "DistanceGeometry/RefinementDebugData.h"
#include "Log.h"
#include "boost_outcome/outcome.hpp"

/*! @file
 *
 * Declares the central conformation (and -ensemble) generating functions that
 * start the DG procedure.
 */

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

namespace predicates {

bool hasZeroPermutationsStereocenters(const Molecule& molecule);

bool hasUnassignedStereocenters(const Molecule& mol);

} // namespace predicates

namespace detail {

AngstromWrapper convertToAngstromWrapper(
  const Eigen::VectorXd& vectorizedPositions
);

AngstromWrapper convertToAngstromWrapper(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
);

/*!
 * A logging, not throwing otherwise identical implementation of
 * runDistanceGeometry, that returns detailed intermediate data from a
 * refinement, while runDistanceGeometry returns only the final result.
 *
 * @note Contained PositionCollections are in Angstrom length units
 */
std::list<RefinementData> debugDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const Partiality& metrizationOption = Partiality::FourAtom,
  const bool& useYInversionTrick = true
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
> runDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const Partiality& metrizationOption = Partiality::FourAtom,
  const bool& useYInversionTrick = true
);

} // namespace detail


//! Intermediate conformational data about a Molecule given by a spatial model
struct MoleculeDGInformation {
  DistanceBoundsMatrix bounds;
  std::vector<ChiralityConstraint> chiralityConstraints;
};

//! Collects intermediate conformational data about a Molecule using a spatial model
MoleculeDGInformation gatherDGInformation(const Molecule& molecule);

/*! Generate a conformational ensemble of a Molecule
 *
 * Returns a result type which may or may not contain a vector of
 * PositionCollections. The result type is much like an optional, except that
 * in the error case it carries data about the error in order to help diagnose
 * possible mistakes made in the molecular graph specification.
 */
outcome::result<
  std::vector<Delib::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  const unsigned numStructures
);

/*! Generate a 3D structure of a Molecule
 *
 * Returns a result type which may or may not contain a PositionCollection. The
 * result type is much like an optional, except that in the error case it
 * carries data about the error in order to help diagnose possible mistakes
 * made in the molecular graph specification.
 */
outcome::result<Delib::PositionCollection> generateConformation(const Molecule& molecule);

} // namespace DistanceGeometry

} // namespace molassembler

#endif
