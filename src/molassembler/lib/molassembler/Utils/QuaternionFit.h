/**
 * @file
 * @brief Slightly modified (to reduce dependencies) file stolen from UtilsOS
 * @note When re-basing from Delib to UtilsOS, this can be removed
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *   See LICENSE.txt for details.
 */
#ifndef UTILS_QUATERNIONFIT_H
#define UTILS_QUATERNIONFIT_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @class QuaternionFit QuaternionFit.h
 * @brief Compares two sets of positions and calculates translation and rotation to go from one to the other.
 *
 * References:
 *  - https://arxiv.org/pdf/physics/0506177.pdf (arXiv:physics/0506177v3)
 *  - https://dx.doi.org/10.1002/jcc.20110
 */
class QuaternionFit {
 public:
  /**
   * @brief Construct a new Quaternion Fit object.
   * @param refMat The positions of the reference data.
   * @param fitMat The positions of the data to be fitted.
   * @param weights The individual weights for each point in the set of data.
   * @param improperRotationIsAllowed If true allows the algorithm to invert the geometry to generate a better fit.
   */
  template<typename DerivedA, typename DerivedB>
  QuaternionFit(const Eigen::DenseBase<DerivedA>& refMat, const Eigen::DenseBase<DerivedB>& fitMat,
                const Eigen::VectorXd& weights, bool improperRotationIsAllowed = false)
    : weights_(weights), refMat_(refMat), fitMat_(fitMat), improperRotationIsAllowed_(improperRotationIsAllowed) {
    assert(refMat_.cols() == 3 && fitMat_.cols() == 3);
    assert(refMat_.rows() == fitMat_.rows());
    assert(refMat_.rows() == weights_.size());
    align();
  };
  /**
   * @brief Construct a new Quaternion Fit object
   * @param Rref The positions of the reference data.
   * @param Rfit The positions of the data to be fitted.
   * @param improperRotationIsAllowed If true allows the algorithm to invert the geometry to generate a better fit.
   */
  template<typename DerivedA, typename DerivedB>
  QuaternionFit(const Eigen::DenseBase<DerivedA>& refMat, const Eigen::DenseBase<DerivedB>& fitMat,
                bool improperRotationIsAllowed = false)
    : QuaternionFit(refMat, fitMat, Eigen::VectorXd::Ones(fitMat.rows()), improperRotationIsAllowed){};
  /**
   * @brief Getter for the reverse of the applied rotation.
   * @return Eigen::Matrix3d The rotation from the reference structure to the second structure before fitting.
   */
  Eigen::Matrix3d getRotationMatrix() const;
  /**
   * @brief Getter for the reverse of the applied translation.
   * @return Eigen::Vector3d The translation from the reference structure to the second structure before fitting.
   */
  Eigen::Vector3d getTransVector() const;
  /**
   * @brief Getter for the fitted data as Eigen3 object.
   * @return Eigen::Matrix3d The fitted data.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> getFittedData() const;
  /**
   * @brief Getter for the RMSD not using any weights that might be stored.
   * @return double The RMSD.
   */
  double getRotRMSD() const;
  /**
   * @brief Getter for the RMSD due to differences in rotation only.
   * @return double The RMSD.
   */
  double getRMSD() const;
  /**
   * @brief Getter for the RMSD using the given weights.
   * @param weights The individaul weights for each point in the set of data.
   * @return double The RMSD.
   */
  double getWeightedRMSD(const Eigen::VectorXd& weights) const;
  /**
   * @brief Getter for the RMSD using the internal weights given/implied in the constructor.
   * @return double The RMSD.
   */
  double getWeightedRMSD() const;

 private:
  void align();

  Eigen::VectorXd weights_;
  Eigen::MatrixX3d refMat_;
  Eigen::MatrixX3d fitMat_;
  Eigen::Vector3d refCenter_;
  Eigen::Vector3d fitCenter_;
  Eigen::Matrix3d rotMat_;
  Eigen::MatrixX3d fittedMat_;
  double maxEigenvalue_;
  bool improperRotationIsAllowed_;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_QUATERNIONFIT_H
