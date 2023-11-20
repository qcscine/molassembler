/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Basic Newton-Raphson minimizer
 */

#ifndef INCLUDE_TEMPLE_NEWTON_RAPHSON_H
#define INCLUDE_TEMPLE_NEWTON_RAPHSON_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/NumericalDiff>
#include <iostream>

namespace Scine {
namespace Molassembler {
namespace Temple {

/**
 * @brief A very basic newton-raphson minimizer
 *
 * No hessian adjustments, no line searches, no trust region safeties, etc.
 * If this minimizes or just approaches saddle points or maxima is
 * uncertain.
 *
 * @tparam FloatType Floating point type of the objective function
 */
template<typename FloatType = double>
struct NewtonRaphson {
  using MatrixType = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;

  //! Type returned from an optimization
  struct OptimizationReturnType {
    //! Number of iterations
    unsigned iterations;
    //! Final function value
    FloatType value;
    //! Final gradient
    VectorType gradient;
  };

  template<
    typename UpdateFunction,
    typename Checker
  > OptimizationReturnType minimize(
    Eigen::Ref<VectorType> parameters,
    UpdateFunction&& function,
    Checker&& check
  ) {
    const unsigned P = parameters.size();
    FloatType value;
    VectorType gradients = Eigen::VectorXd::Zero(P);
    MatrixType hessian (P, P);
    hessian.setZero();

    VectorType numericalGradient(P);
    {
      constexpr FloatType h = 1e-4;
      FloatType a;
      FloatType b;
      VectorType diffParameters = parameters;
      for(unsigned i = 0; i < P; ++i) {
        diffParameters(i) += h;
        function(diffParameters, b, gradients, hessian);
        diffParameters(i) -= 2 * h;
        function(diffParameters, a, gradients, hessian);
        diffParameters(i) = parameters(i);
        numericalGradient(i) = (b - a) / (2 * h);
      }
    }

    function(parameters, value, gradients, hessian);

    if(numericalGradient.norm() > 1e-4) {
      if(!numericalGradient.isApprox(gradients, 1e-8)) {
        std::cout << "MISMATCH: Numerical gradient is " << numericalGradient.transpose() << "\nfunction gradient is : " << gradients.transpose() << "\n";
      } else {
        std::cout << "match: norm is " << gradients.norm() << " and diff norm is " << (numericalGradient - gradients).norm() << "\n";
      }
    }

    MatrixType numericalHessian(P, P);
    {
      constexpr FloatType h = 1e-4;
      FloatType a;
      FloatType b;
      FloatType c;
      FloatType d;
      VectorType diffParameters = parameters;
      for(unsigned i = 0; i < P; ++i) {
        // Diagonal second derivative from forward and backward derivative
        function(diffParameters, b, gradients, hessian);
        diffParameters(i) += h;
        function(diffParameters, c, gradients, hessian);
        diffParameters(i) -= 2 * h;
        function(diffParameters, a, gradients, hessian);
        diffParameters(i) = parameters(i);
        numericalHessian(i, i) = (a - 2 * b + c) / (h * h);

        // Cross terms from two central derivatives
        for(unsigned j = i + 1; j < P; ++j) {
          // (i) + h, (j) + h -> a
          diffParameters(i) += h;
          diffParameters(j) += h;
          function(diffParameters, a, gradients, hessian);
          // (i) + h, (j) - h - > b
          diffParameters(j) -= 2 * h;
          function(diffParameters, b, gradients, hessian);
          // (i) - h, (j) - h -> d
          diffParameters(i) -= 2 * h;
          function(diffParameters, d, gradients, hessian);
          // (i) - h, (j) + h -> c
          diffParameters(j) += 2 * h;
          function(diffParameters, c, gradients, hessian);

          numericalHessian(i, j) = (a - b - c + d) / (4 * h * h);
          numericalHessian(j, i) = numericalHessian(i, j);
        }
      }
    }

    function(parameters, value, gradients, hessian);

    if(!numericalHessian.isApprox(hessian), 1e-2) {
      std::cout << "MISMATCH: Hessian is\n" << hessian << "\nNumerical hessian is\n" << numericalHessian << "\n";
    } else {
      std::cout << "Hessian matches\n";
    }

    unsigned iteration = 0;
    while(
      check.shouldContinue(
        iteration,
        std::as_const(value),
        std::as_const(gradients)
      )
    ) {
      // Solve HΔx = g, then apply x_(n+1) = x_n - Δx
      Eigen::JacobiSVD<MatrixType> svd(
        hessian,
        Eigen::ComputeFullV | Eigen::ComputeFullU
      );
      svd.setThreshold(svdThreshold);
      VectorType steps = svd.solve(gradients);

      // Take a step
      FloatType rms = std::sqrt(steps.squaredNorm() / P);
      if (rms > trustRadius) {
        steps *= trustRadius / rms;
      }
      parameters -= steps;

      function(parameters, value, gradients, hessian);
      ++iteration;
    }

    return {
      iteration,
      value,
      std::move(gradients)
    };
  }

  //! The SVD threshold for the decomposition of the Hessian
  FloatType svdThreshold = 1.0e-12;
  //! The maximum RMS of a taken step
  FloatType trustRadius = 0.5;
};

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
