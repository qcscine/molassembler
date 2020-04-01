/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Shapes/InertialMoments.h"

#include "molassembler/Shapes/CoordinateSystemTransformation.h"
#include "molassembler/Shapes/ContinuousMeasures.h"
#include <Eigen/Eigenvalues>

#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/Adaptors/Iota.h"

namespace Scine {
namespace shapes {
namespace detail {
//! Determine degeneracy of intertial moments
unsigned degeneracy(const Eigen::Vector3d& inertialMoments) {
  constexpr double degeneracyEpsilon = 0.05;
  unsigned mdeg = 0;
  if(
    std::fabs(
      (inertialMoments(2) - inertialMoments(1)) / inertialMoments(2)
    ) <= degeneracyEpsilon
  ) {
    mdeg += 1;
  }
  if(
    std::fabs(
      (inertialMoments(1) - inertialMoments(0)) / inertialMoments(2)
    ) <= degeneracyEpsilon
  ) {
    mdeg += 2;
  }

  return 1 + (mdeg + 1) / 2;
}

} // namespace detail

InertialMoments principalInertialMoments(
  const InertialPositionsType& normalizedPositions
) {
  Eigen::Matrix3d inertialMatrix = Eigen::Matrix3d::Zero(3, 3);
  const unsigned N = normalizedPositions.cols();

  for(unsigned i = 0; i < N; ++i) {
    const auto& vec = normalizedPositions.col(i);
    inertialMatrix(0, 0) += vec.y() * vec.y() + vec.z() * vec.z();
    inertialMatrix(1, 1) += vec.x() * vec.x() + vec.z() * vec.z();
    inertialMatrix(2, 2) += vec.x() * vec.x() + vec.y() * vec.y();
    inertialMatrix(1, 0) -= vec.x() * vec.y(); // xy
    inertialMatrix(2, 0) -= vec.x() * vec.z(); // xz
    inertialMatrix(2, 1) -= vec.y() * vec.z(); // yz
  }

  // Decompose the inertial matrix to get principal axes and inertial moments
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> decomposition(inertialMatrix);

  InertialMoments result;
  result.moments = decomposition.eigenvalues();
  result.axes = decomposition.eigenvectors();
  return result;
}

Top standardizeTop(Eigen::Ref<InertialPositionsType> normalizedPositions) {
  const unsigned N = normalizedPositions.cols();
  assert(N > 1);

  InertialMoments moments = principalInertialMoments(normalizedPositions);

  const unsigned degeneracy = detail::degeneracy(moments.moments);

  auto rotateEverything = [&](const CoordinateSystem& sourceSystem) {
    const CoordinateSystem defaultCoordinateSystem {};
    assert(sourceSystem.isRightHanded());
    const auto R = rotationMatrix(sourceSystem, defaultCoordinateSystem);
    // Rotate coordinates
    normalizedPositions = R * normalizedPositions;
    // Rotate inertial moment axes
    moments.axes = R * moments.axes;
  };

  if(moments.moments(0) < 0.1 && degeneracy == 2) {
    // The top is linear: If IA << IB = IC and IA ~ 0. We rotate IA to z
    const CoordinateSystem inertialMomentSystem {
      moments.axes.col(1),
      moments.axes.col(2)
    };
    rotateEverything(inertialMomentSystem);
    assert(moments.axes.col(0).cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));
    return Top::Line;
  }

  if(degeneracy == 1) {
    /* The top is asymmetric. Rotate the axis with the
     * highest moment of inertia to coincide with z, and the one with second most
     * to coincide with x.
     *
     * To better define orientation, we could look for Cn axes. This is done in
     * another function. No need to burden this function with that here.
     */
    CoordinateSystem inertialMomentSystem {
      moments.axes.col(1), // second highest becomes x
      moments.axes.col(2).cross(moments.axes.col(1)) // y = z.cross(x)
    };
    assert(inertialMomentSystem.z.isApprox(moments.axes.col(2), 1e-10));
    rotateEverything(inertialMomentSystem);
    // Make sure rotation went as intended
    assert(moments.axes.col(2).cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));
    assert(moments.axes.col(1).cwiseAbs().isApprox(Eigen::Vector3d::UnitX(), 1e-10));

    return Top::Asymmetric;
  }

  if(degeneracy == 2) {
    /* The top is symmetric. The subsets are:
     * - Oblate (disc): IA = IB < IC
     * - Prolate (rugby football): IA < IB = IC
     *
     * We rotate the unique axis to coincide with z (it's probably the site of
     * the highest-order Cn or Sn, and one of the degenerate axes to coincide
     * with x. There could be a C2 on x.
     *
     * This is most likely rare and should occur only for largely undistorted
     * structures. Perhaps we can flowchart point groups here?
     */
    // Calculate Ray's asymmetry parameter
    const double A = 1 / moments.moments(0);
    const double B = 1 / moments.moments(1);
    const double C = 1 / moments.moments(2);
    const double kappa = (2 * B - A - C) / (A - C);
    assert(-1 <= kappa && kappa <= 1);
    if(kappa < 0) {
      // Prolate top. IA is unique
      CoordinateSystem inertialMomentSystem {
        moments.axes.col(1),
        moments.axes.col(2)
      };
      rotateEverything(inertialMomentSystem);
      assert(moments.axes.col(0).cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));
      return Top::Prolate;
    }

    // Oblate top. IC is unique
    CoordinateSystem inertialMomentSystem {
      moments.axes.col(0),
      moments.axes.col(1)
    };
    rotateEverything(inertialMomentSystem);
    assert(moments.axes.col(2).cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));
    return Top::Oblate;
  }

  assert(degeneracy == 3);
  /* The top is spherical (IA = IB = IC).
   *
   * Note that there is no reason to rotate anything on the basis of the axes
   * from the inertial moment analysis since any choice of axes gives the
   * spherical symmetry. We can't use those to rotate the system.
   *
   * Rotate an arbitrary position to +z instead (good for Td
   * and Oh octahedral, less so for Oh cubic and Ih, which should be less
   * common)
   */
  unsigned selectedIndex = 0;
  for(; selectedIndex < N; ++selectedIndex) {
    /* As long as the position isn't close to the centroid and it's not exactly
     * the -z vector, we can rotate it
     */
    if(
      normalizedPositions.col(selectedIndex).norm() > 0.2
      && !normalizedPositions.col(selectedIndex).normalized().isApprox(
        -Eigen::Vector3d::UnitZ(),
        1e-10
      )
    ) {
      break;
    }
  }
  assert(selectedIndex != N);

  // Determine axis of rotation as sum of z and position coordinate
  const Eigen::Vector3d rotationAxis = (
    normalizedPositions.col(selectedIndex).normalized()
    + Eigen::Vector3d::UnitZ()
  ).normalized();
  const Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(M_PI, rotationAxis).toRotationMatrix();

  // Rotate all coordinates
  for(unsigned i = 0; i < N; ++i) {
    normalizedPositions.col(i) = rotationMatrix * normalizedPositions.col(i);
  }
  // Check that everything went as planned
  assert(normalizedPositions.col(selectedIndex).normalized().cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));

  return Top::Spherical;
}

unsigned reorientAsymmetricTop(Eigen::Ref<InertialPositionsType> normalizedPositions) {
  const unsigned P = normalizedPositions.cols();
  const auto& axes = Eigen::Matrix3d::Identity();

  struct AxisBest {
    unsigned order = 1;
    double csm = 1; // This functions much like a detection threshold below
    unsigned axisIndex;

    AxisBest(unsigned index) : axisIndex(index) {}

    bool operator < (const AxisBest& other) const {
      return order > other.order;
    }
  };

  auto orderedAxisBest = temple::sort(
    temple::map(
      temple::iota<unsigned>(3),
      [&](const unsigned axisIndex) -> AxisBest {
        const Eigen::Vector3d axis = axes.col(axisIndex);
        AxisBest best {axisIndex};
        for(unsigned n = 2; n <= P; ++n) {
          const double axisCSM = continuous::fixed::element(normalizedPositions, elements::Rotation::Cn(axis, n));
          if(axisCSM < best.csm) {
            best.order = n;
            best.csm = axisCSM;
          }
        }
        return best;
      }
    )
  );

  if(orderedAxisBest.front().order > 1) {
    /* Only mess with the coordinate frame if any sort of axis was found.
     * We want the second-highest order axis on x, highest order axis along z,
     * doesn't really matter if +z or -z
     */
    const CoordinateSystem highestOrderSystem {
      axes.col(orderedAxisBest.at(1).axisIndex),
      axes.col(orderedAxisBest.back().axisIndex)
    };

    normalizedPositions = rotationMatrix(highestOrderSystem, {}) * normalizedPositions;
  }

  return orderedAxisBest.front().order;
}

} // namespace shapes
} // namespace Scine
