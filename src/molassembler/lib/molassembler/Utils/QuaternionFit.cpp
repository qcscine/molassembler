/**
 * @file
 * @brief Slightly modified (to reduce dependencies) file stolen from UtilsOS
 * @note When re-basing from Delib to UtilsOS, this can be removed
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *   See LICENSE.txt for details.
 */
#include "molassembler/Utils/QuaternionFit.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <utility>

using namespace Eigen;

namespace Scine {
namespace Utils {

void QuaternionFit::align() {
  /*
   * Translation into origin
   */
  // multiply positions with weights colwise: mat.array().colwise()*weights.array()
  // sum into vector and divide: (...).sum() / weights.sum()
  refCenter_ = (refMat_.array().colwise() * weights_.array()).colwise().sum() / weights_.sum();
  fitCenter_ = (fitMat_.array().colwise() * weights_.array()).colwise().sum() / weights_.sum();
  fittedMat_ = fitMat_.rowwise() - fitCenter_.transpose();

  /*
   * Rotation
   */
  Eigen::Matrix4d b = Eigen::Matrix4d::Zero();
  // generate decomposable matrix per atom and add them
  for (int i = 0; i < fitMat_.rows(); i++) {
    Eigen::Vector3d fitPos = fittedMat_.row(i);
    Eigen::Vector3d refPos = refMat_.row(i) - refCenter_.transpose();
    Eigen::Matrix4d a = Eigen::Matrix4d::Zero();
    a.block(0, 1, 1, 3) = (fitPos - refPos).transpose();
    a.block(1, 0, 3, 1) = refPos - fitPos;
    a.block(1, 1, 3, 3) = Eigen::Matrix3d::Identity().rowwise().cross(refPos + fitPos);
    b += a.transpose() * a * weights_[i];
  }

  // Decompose b
  SelfAdjointEigenSolver<Matrix4d> eigensolver(b);

  // Apply rotation
  //   If the eigenvalue of the last eigenvector is larger (absolute) than that of the first one
  //   it is beneficial to allow a rotation including inversion, as this will lead to the better
  //   fit.
  if (improperRotationIsAllowed_ && fabs(eigensolver.eigenvalues()[0]) < fabs(eigensolver.eigenvalues()[3])) {
    const Eigen::Vector4d& q = eigensolver.eigenvectors().col(3);
    rotMat_ = -Quaterniond(q[0], q[1], q[2], q[3]).toRotationMatrix();
    fittedMat_ = (rotMat_ * fittedMat_.transpose()).transpose();
    maxEigenvalue_ = eigensolver.eigenvalues()[3];
  }
  else {
    const Eigen::Vector4d& q = eigensolver.eigenvectors().col(0);
    rotMat_ = Quaterniond(q[0], q[1], q[2], q[3]).toRotationMatrix();
    fittedMat_ = (rotMat_ * fittedMat_.transpose()).transpose();
    maxEigenvalue_ = eigensolver.eigenvalues()[0];
  }

  /*
   * Translation onto reference center
   */
  fittedMat_ = fittedMat_.rowwise() + refCenter_.transpose();
}

double QuaternionFit::getRMSD() const {
  // collect squared norms: ((refMat_-fittedMat_).rowwise().squaredNorm()
  // sum, divide and sqrt: std::sqrt( (...).sum() / refMat_.rows() )
  return std::sqrt((refMat_ - fittedMat_).rowwise().squaredNorm().sum() / refMat_.rows());
}

double QuaternionFit::getWeightedRMSD(const Eigen::VectorXd& weights) const {
  // collect sqared norms: ((refMat_-fittedMat_).rowwise().squaredNorm()
  // multiply with corresponding weight: (...).array()*weights.array()
  // sum, divide and sqrt: std::sqrt( (...).sum() / refMat_.rows() )
  return std::sqrt(((refMat_ - fittedMat_).rowwise().squaredNorm().array() * weights.array()).sum() / refMat_.rows());
}

double QuaternionFit::getWeightedRMSD() const {
  return std::move(getWeightedRMSD(weights_));
}

Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> QuaternionFit::getFittedData() const {
  return fittedMat_;
}

Eigen::Matrix3d QuaternionFit::getRotationMatrix() const {
  return rotMat_.transpose();
}

Eigen::Vector3d QuaternionFit::getTransVector() const {
  return fitCenter_ - refCenter_;
}

double QuaternionFit::getRotRMSD() const {
  double rotRMSD = (fitMat_.rowwise() - fitCenter_.transpose()).rowwise().squaredNorm().sum();
  rotRMSD += (refMat_.rowwise() - refCenter_.transpose()).rowwise().squaredNorm().sum();
  rotRMSD -= 2.0 * abs(maxEigenvalue_);
  if (rotRMSD > 0)
    return sqrt(rotRMSD / refMat_.rows());
  return 0.0;
}

} /* namespace Utils */
} /* namespace Scine */
