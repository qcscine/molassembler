/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Data types common to optimizers
 */

#ifndef INCLUDE_TEMPLE_OPTIMIZATION_COMMON_H
#define INCLUDE_TEMPLE_OPTIMIZATION_COMMON_H

#include <Eigen/Core>

namespace Scine {
namespace Temple {

//! @brief Functionality common to multiple optimizers
namespace Optimization {

template<typename FloatType>
struct FloatUpdateBuffer {
  static_assert(
    std::is_floating_point<FloatType>::value,
    "This struct is not intended for anything else but floats"
  );

  FloatType current;
  FloatType proposed;

  //! Returns signed difference: proposed - current
  FloatType absDelta() const {
    return std::fabs(proposed - current);
  }

  void propagate() { current = proposed; }
};

template<typename EigenType>
struct EigenUpdateBuffer {
  EigenType current;
  EigenType proposed;

  typename EigenType::Scalar deltaNorm() const {
    return (proposed - current).norm();
  }

  void propagate() {
    /* Swapping instead of moving proposed into current is three moves instead
     * of one, but avoids allocation of new memory (after moving from proposed,
     * we need to resize proposed to have the same size). Assuming moves are
     * cheap (trade ownership of allocated memory) and allocation could be
     * costly, this might be cheaper.
     */
    std::swap(current, proposed);
  }
};

template<typename FloatType, typename VectorType, typename UpdateFunction>
struct NegateFunction {
  UpdateFunction function;

  NegateFunction(UpdateFunction&& fn) : function(fn) {}

  void operator() (
    const VectorType& parameters,
    FloatType& value,
    Eigen::Ref<VectorType> gradient
  ) {
    function(parameters, value, gradient);
    value *= -1;
    gradient *= -1;
  }
};

template<typename VectorType, typename UpdateFunction>
auto negateFunction(UpdateFunction&& fn) {
  using FloatType = typename VectorType::Scalar;
  return NegateFunction<FloatType, VectorType, UpdateFunction>(std::forward<UpdateFunction>(fn));
}

template<typename Function, typename FloatType>
auto numericalGradient(
  Function&& function,
  const Eigen::Matrix<FloatType, Eigen::Dynamic, 1>& parameters
) {
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;
  const unsigned P = parameters.size();
  VectorType numericalGradient(P);
  VectorType diffParameters = parameters;

  // Central differences for all individual parameters
  constexpr FloatType h = 1e-4;
  FloatType a, b;
  for(unsigned i = 0; i < P; ++i) {
    diffParameters(i) += h;
    b = function(diffParameters);
    diffParameters(i) -= 2 * h;
    a = function(diffParameters);
    diffParameters(i) = parameters(i);
    numericalGradient(i) = (b - a) / (2 * h);
  }

  return numericalGradient;
}

template<typename Function, typename FloatType>
auto numericalHessian(
  Function&& function,
  const Eigen::Matrix<FloatType, Eigen::Dynamic, 1>& parameters
) {
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;
  using MatrixType = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>;
  const unsigned P = parameters.size();
  MatrixType numericalHessian(P, P);
  VectorType diffParameters = parameters;

  // Form the hessian using finite differences
  constexpr FloatType h = 1e-4;
  FloatType a, b, c, d;
  for(unsigned i = 0; i < P; ++i) {
    // Diagonal second derivative from forward and backward derivative
    b = function(diffParameters);
    diffParameters(i) += h;
    c = function(diffParameters);
    diffParameters(i) -= 2 * h;
    a = function(diffParameters);
    diffParameters(i) = parameters(i);
    numericalHessian(i, i) = (a - 2 * b + c) / (h * h);

    // Cross terms from two central derivatives
    for(unsigned j = i + 1; j < P; ++j) {
      // (i) + h, (j) + h -> a
      diffParameters(i) += h;
      diffParameters(j) += h;
      a = function(diffParameters);
      // (i) + h, (j) - h - > b
      diffParameters(j) -= 2 * h;
      b = function(diffParameters);
      // (i) - h, (j) - h -> d
      diffParameters(i) -= 2 * h;
      d = function(diffParameters);
      // (i) - h, (j) + h -> c
      diffParameters(j) += 2 * h;
      c = function(diffParameters);

      numericalHessian(i, j) = (a - b - c + d) / (4 * h * h);
      numericalHessian(j, i) = numericalHessian(i, j);
    }
  }

  return numericalHessian;
}

} // namespace Optimization
} // namespace Temple
} // namespace Scine

#endif
