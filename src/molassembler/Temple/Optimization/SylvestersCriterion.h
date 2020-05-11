/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Sylvester's criterion for positive definiteness and positive
 *   semi-definiteness
 */

#ifndef INCLUDE_TEMPLE_SYLVESTERS_CRITERION_H
#define INCLUDE_TEMPLE_SYLVESTERS_CRITERION_H

#include <Eigen/Core>

namespace Scine {
namespace Molassembler {
namespace Temple {

/**
 * @brief Determine whether a matrix is positive definite
 *
 * @pre @p matrix must be hermitian
 *
 * @return
 */
template<typename Derived>
bool positiveDefinite(const Eigen::MatrixBase<Derived>& matrix) {
  using Scalar = typename Eigen::MatrixBase<Derived>::Scalar;
  assert(matrix.rows() == matrix.cols());

  const unsigned N = matrix.rows();

  // 1st leading principal minor
  if(matrix(0, 0) <= Scalar {0.0}) {
    return false;
  }

  if(N == 1) {
    return true;
  }

  /* Test leading principal minors of sizes 2-4 using fixed-size matrix
   * subexpressions, those should be significantly faster
   */
  // 2nd leading principal minor
  if(matrix.template block<2, 2>(0, 0).determinant() <= Scalar {0.0}) {
    return false;
  }

  if(N == 2) {
    return true;
  }

  // 3rd leading principal minor
  if(matrix.template block<3, 3>(0, 0).determinant() <= Scalar {0.0}) {
    return false;
  }

  if(N == 3) {
    return true;
  }

  // 4th leading principal minor
  if(matrix.template block<4, 4>(0, 0).determinant() <= Scalar {0.0}) {
    return false;
  }

  if(N == 4) {
    return true;
  }

  /* There is limited speed advantage to fixed-size matrix sub-expressions from
   * here, so we can now loop over the rest:
   */
  for(unsigned I = 5; I <= N; ++I) {
    if(matrix.block(0, 0, I, I).determinant() <= Scalar {0.0}) {
      return false;
    }
  }

  return true;
}

/**
 * @brief Determine whether a matrix is positive semidefinite
 *
 * @pre @p matrix must be hermitian
 *
 * @return
 */
template<typename Derived>
bool positiveSemidefinite(const Eigen::MatrixBase<Derived>& matrix) {
  using Scalar = typename Eigen::MatrixBase<Derived>::Scalar;
  assert(matrix.rows() == matrix.cols());

  const unsigned N = matrix.rows();

  // Principal minors of size I = 1
  if((matrix.array() < Scalar {0.0}).any()) {
    return false;
  }

  if(N == 1) {
    return true;
  }

  // Principal minors of size I = 2
  for(unsigned i = 0; i < N - 1; ++i) {
    for(unsigned j = 0; j < N - 1; ++j) {
      if(matrix.template block<2, 2>(i, j).determinant() < Scalar {0.0}) {
        return false;
      }
    }
  }

  if(N == 2) {
    return true;
  }

  // Principal minors of size I = 3
  for(unsigned i = 0; i < N - 2; ++i) {
    for(unsigned j = 0; j < N - 2; ++j) {
      if(matrix.template block<3, 3>(i, j).determinant() < Scalar {0.0}) {
        return false;
      }
    }
  }

  if(N == 3) {
    return true;
  }

  // Principal minors of size I = 4
  for(unsigned i = 0; i < N - 3; ++i) {
    for(unsigned j = 0; j < N - 3; ++j) {
      if(matrix.template block<4, 4>(i, j).determinant() < Scalar {0.0}) {
        return false;
      }
    }
  }

  if(N == 4) {
    return true;
  }

  // Remaining principal minors
  for(unsigned I = 5; I <= N; ++I) {
    for(unsigned i = 0; i < N - I; ++i) {
      for(unsigned j = 0; j < N - I; ++j) {
        if(matrix.block(i, j, I, I).determinant() < Scalar {0.0}) {
          return false;
        }
      }
    }
  }

  return true;
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
