/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Distance Geometry conformation generating procedures
 *
 * Declares the central conformation (and -ensemble) generating functions that
 * start the DG procedure.
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_CONFORMER_GENERATION_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_CONFORMER_GENERATION_H

#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/Log.h"

namespace Scine {

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

namespace detail {

/*! @brief Collects four-dimensional linear positions into three-dimensional matrix
 *
 * Vectorized positions are x0, y0, z0, w0, x1, y1, z1, w1, ...
 * Result of gather is Nx3 matrix
 *
 * @complexity{@math{\Theta(N)}}
 */
Eigen::MatrixXd gather(const Eigen::VectorXd& vectorizedPositions);

/**
 * @brief Converts Nx3 matrix into AngstromWrapper type
 *
 * @complexity{@math{\Theta(N)}}
 */
AngstromWrapper convertToAngstromWrapper(const Eigen::MatrixXd& positions);

/**
 * @brief Rotates and translates the generated coordinates against the fixed
 *   positions
 *
 * @complexity{A single quaternion fit, possibly linear in the number of fixed
 * positions or linear in the number of atoms}
 */
Eigen::MatrixXd fitAndSetFixedPositions(
  const Eigen::MatrixXd& positions,
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
 *
 * @complexity{At least linear in the number of unassigned stereopermutators
 * multiplied by the number of atoms.}
 */
Molecule narrow(Molecule molecule, random::Engine& engine);

} // namespace detail

//! Intermediate conformational data about a Molecule given by a spatial model
struct MoleculeDGInformation {
  SpatialModel::BoundsMatrix bounds;
  std::vector<ChiralConstraint> chiralConstraints;
  std::vector<DihedralConstraint> dihedralConstraints;
};

/*! @brief Collects intermediate conformational data about a Molecule using a spatial model
 *
 * @complexity{At least @math{O(P_2 + P_3 + P_4)} where @math{P_i} is the
 * number of distinct paths of length @math{i} in the graph. That should
 * scale at least linearly in the number of vertices.}
 */
MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const Configuration& configuration
);

//! @brief Distance Geometry refinement
outcome::result<AngstromWrapper> refine(
  Eigen::MatrixXd embeddedPositions,
  const DistanceBoundsMatrix& distanceBounds,
  const Configuration& configuration,
  const std::shared_ptr<MoleculeDGInformation>& DgDataPtr
);

// @brief Individual conformer generation routine
outcome::result<AngstromWrapper> generateConformer(
  const Molecule& molecule,
  const Configuration& configuration,
  std::shared_ptr<MoleculeDGInformation>& DgDataPtr,
  bool regenerateDGDataEachStep,
  random::Engine& engine
);

/** @brief Main and parallel implementation of Distance Geometry. Generates an
 *   ensemble of 3D structures of a given Molecule
 *
 * @complexity{Roughly @math{O(C \cdot N^3) where @math{C} is the number of
 * conformers and @math{N} is the number of atoms in @p molecule}
 *
 * @see generateEnsemble
 */
std::vector<
  outcome::result<AngstromWrapper>
> run(
  const Molecule& molecule,
  unsigned numConformers,
  const Configuration& configuration,
  const boost::optional<unsigned> seedOption
);

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine

#endif
