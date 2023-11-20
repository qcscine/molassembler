/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "CoordinateSystemTransformation.h"

namespace Scine {
namespace Molassembler {
namespace Shapes {
namespace Detail {

inline double angle(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
  return std::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

inline double signedAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& n) {
  return std::copysign(
    angle(a, b),
    n.dot(a.cross(b))
  );
}

} // namespace Detail

CoordinateSystem::CoordinateSystem() : x(Eigen::Vector3d::UnitX()), y(Eigen::Vector3d::UnitY()), z(Eigen::Vector3d::UnitZ()) {}
CoordinateSystem::CoordinateSystem(const Eigen::Vector3d& a, const Eigen::Vector3d& b) : x(a.normalized()), y(b.normalized()), z(x.cross(y).normalized()) {
  // a and b need to be at right angles
  assert(std::fabs(Detail::angle(a, b) - M_PI / 2) < 1e-4);
}

CoordinateSystem CoordinateSystem::random() {
  const Eigen::Vector3d a = Eigen::Vector3d::Random().normalized();
  return {
    a,
    Eigen::Vector3d::Random().cross(a).normalized()
  };
}

bool CoordinateSystem::isRightHanded() const {
  return x.normalized().cross(y.normalized()).isApprox(z.normalized(), 1e-10);
}

bool isRotationMatrix(const Eigen::Matrix3d& R) {
  return R.transpose().isApprox(R.inverse(), 1e-10);
}

Eigen::Matrix3d rotationMatrix(const CoordinateSystem& a, const CoordinateSystem& b) {
  assert(a.isRightHanded());
  assert(b.isRightHanded());

  /* If z is z', then z x z' is a null vector, with which we cannot rotate
   * anything. Catch the case that z matches already:
   */
  if(a.z.isApprox(b.z, 1e-8)) {
    Eigen::Matrix3d R = Eigen::AngleAxisd(
      Detail::signedAngle(a.x, b.x, a.z),
      a.z.normalized()
    ).toRotationMatrix();
    assert(isRotationMatrix(R));
    return R;
  }

  /* If z is -z', then z x z' is still a null vector with which we cannot
   * rotate anything. Catch that case too:
   */
  if(a.z.isApprox(-b.z, 1e-8)) {
    Eigen::Matrix3d R = -Eigen::AngleAxisd(
      -Detail::signedAngle(a.x, b.x, a.z),
      a.z.normalized()
    ).toRotationMatrix();
    assert(isRotationMatrix(R));
    return R;
  }

  const Eigen::Vector3d N = a.z.cross(b.z);
  /* Angle definitions:
   * - alpha: x and N viewed from z
   * - beta: z and z' viewed from N
   * - gamma: N and x' viewed from z'
   */
  const double alpha = Detail::signedAngle(a.x, N, a.z);
  const double beta = Detail::signedAngle(a.z, b.z, N);
  const double gamma = Detail::signedAngle(N, b.x, b.z);

  assert(!std::isnan(alpha) && !std::isnan(beta) && !std::isnan(gamma));

  Eigen::Matrix3d R = (
    Eigen::AngleAxisd(gamma, b.z.normalized())
    * Eigen::AngleAxisd(beta, N.normalized())
    * Eigen::AngleAxisd(alpha, a.z.normalized())
  ).toRotationMatrix();
  assert(isRotationMatrix(R));
  return R;
}

} // namespace Shapes
} // namespace Molassembler
} // namespace Scine
