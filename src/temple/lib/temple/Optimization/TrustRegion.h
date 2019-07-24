/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Trust region minimizer
 */

#ifndef INCLUDE_TEMPLE_OPTIMIZATION_TRUST_REGION_NEWTON_H
#define INCLUDE_TEMPLE_OPTIMIZATION_TRUST_REGION_NEWTON_H

#include "temple/STL17.h"
#include "temple/Optimization/Common.h"
#include "temple/Optimization/SylvestersCriterion.h"
#include "boost/math/tools/roots.hpp"
#include <Eigen/Eigenvalues>

#include <iostream>

namespace temple {

namespace detail {

template<typename Derived>
auto resolve(const Eigen::MatrixBase<Derived>& matrix) -> typename Eigen::MatrixBase<Derived>::Scalar {
  assert(matrix.rows() == 1 && matrix.cols() == 1);
  return matrix(0, 0);
}

template<typename FloatType>
struct TROFloatingPointTolerances {};

template<>
struct TROFloatingPointTolerances<double> {
  static constexpr unsigned rootBitAccuracy = 20;
};

template<>
struct TROFloatingPointTolerances<float> {
  static constexpr unsigned rootBitAccuracy = 10;
};

} // namespace detail

/**
 * @brief Trust region Newton-Raphson minimizer with 2D subspace minimization
 *
 * @tparam FloatType floating point type of the objective function
 */
template<typename FloatType = double>
struct TrustRegionOptimizer {
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
  > static OptimizationReturnType minimize(
    Eigen::Ref<VectorType> parameters,
    UpdateFunction&& function,
    Checker&& check
  ) {
    // This is chosen arbitrarily
    constexpr FloatType trustRadiusLimit = 4.0;
    static_assert(trustRadiusLimit > 0, "Trust radius is not positive!");

    // This is a fixed constant of the algorithm
    constexpr FloatType agreementContractionUpperBound = 0.25;
    constexpr FloatType agreementExpansionLowerBound = 0.75;

    // This is chosen arbitrarily within [0, agreementContractionUpperBound)
    constexpr FloatType agreementEta = 0.125;
    static_assert(
      0.0 <= agreementEta && agreementEta < 0.25,
      "Eta not within bounds"
    );

    // Loop variables initialization
    // Choose initial trust radius within (0, trustRadiusLimit)
    FloatType trustRadius = 1.0; // Δ
    FloatType modelAgreement = 0; // ρ
    unsigned iterations = 0;
    StepValues step;
    step.initialize(parameters, function, trustRadius);

    for(
      ;
      check.shouldContinue(
        iterations,
        stl17::as_const(step.values.current),
        stl17::as_const(step.gradients.current)
      );
      ++iterations
    ) {
      modelAgreement = step.modelAgreement();

      /* Adjust trust radius if necessary */
      if(modelAgreement < agreementContractionUpperBound) {
        // Contract the trust radius
        trustRadius = agreementContractionUpperBound * trustRadius;
      } else if(
        modelAgreement > agreementExpansionLowerBound
        && std::fabs(step.direction.norm() - trustRadius) < 1e-5
      ) {
        // Expand the trust radius
        trustRadius = std::min(2 * trustRadius, trustRadiusLimit);
      }

      /* Decide what to do with the step */
      if(modelAgreement > agreementEta) {
        // Accept the step and prepare a new one
        step.accept(function, trustRadius);
      } else {
        // Keep the current parameters, generate new prospective parameters
        step.prepare(function, trustRadius);
      }
    }

    return {
      iterations,
      step.values.current,
      std::move(step.gradients.current)
    };
  }

private:
  /**
   * @brief Stores all parameter, value, gradient and hessian data of two points
   *
   * Stores the parameters, an increment direction, the resulting proposed
   * parameters, and the function evaluation results for both sets of
   * parameters.
   */
  struct StepValues {
    optimization::EigenUpdateBuffer<VectorType> parameters;
    optimization::FloatUpdateBuffer<FloatType> values;
    optimization::EigenUpdateBuffer<VectorType> gradients;
    optimization::EigenUpdateBuffer<MatrixType> hessians;

    VectorType direction;

    template<typename UpdateFunction>
    void initialize(
      Eigen::Ref<VectorType> initialParameters,
      UpdateFunction&& function,
      const FloatType trustRadius
    ) {
      parameters.current = std::move(initialParameters);
      const unsigned P = parameters.current.size();
      parameters.proposed.resize(P);
      gradients.current.resize(P);
      gradients.proposed.resize(P);
      hessians.current.resize(P, P);
      hessians.proposed.resize(P, P);

      function(parameters.current, values.current, gradients.current, hessians.current);
      prepare(function, trustRadius);
    }

    /**
     * @brief Find the minimizer of the trust region subspace
     *
     * @param trustRadius
     * @param g The gradient vector
     * @param B A positive definite matrix approximating the hessian
     *
     * Solves \f$\min_p m(p) = f + g^Tp + \frac{1}{2}p^TBp\f$
     * subject to \f$||p|| \le \Delta\f$ and \f$p \in \textrm{span}[g,B^{-1}g]\f$
     *
     * @return A direction to increment the parameters by leading to a decrease
     *   in the objective function
     */
    VectorType subspaceMinimize(
      const FloatType trustRadius,
      const VectorType& g,
      const MatrixType& B
    ) const {
      using Vector2Type = Eigen::Matrix<FloatType, 2, 1>;
      using Matrix2Type = Eigen::Matrix<FloatType, 2, 2>;

      assert(positiveDefinite(B));
      const MatrixType Binverse = B.inverse();

      const Vector2Type gTilde {g.squaredNorm(), g.dot(Binverse * g)};
      Matrix2Type BTilde;
      BTilde << g.dot(B * g), gTilde(0),
                gTilde(0), gTilde(1);

      const Vector2Type uStar = - BTilde.inverse() * gTilde;

      Matrix2Type BBar;
      BBar << gTilde(0), gTilde(1),
              gTilde(1), (Binverse * g).squaredNorm();
      BBar *= 2;

      const FloatType J = FloatType {0.5} * uStar.transpose() * BBar * uStar;
      const FloatType trustRadiusSquare = trustRadius * trustRadius;

      if(J <= trustRadiusSquare) {
        return uStar(0) * g + uStar(1) * Binverse * g;
      }

      // Now we have to find lambda
      boost::uintmax_t iterations = 1000;
      auto root_result = boost::math::tools::toms748_solve(
        [&](const FloatType lambda) -> FloatType {
          const Vector2Type intermediate = (BTilde + lambda * BBar).inverse() * gTilde;
          return FloatType {0.5} * intermediate.transpose() * BBar * intermediate - trustRadiusSquare;
        },
        std::nextafter(FloatType {0.0}, FloatType {1.0}),
        FloatType {100.0},
        boost::math::tools::eps_tolerance<FloatType> {
          detail::TROFloatingPointTolerances<FloatType>::rootBitAccuracy
        },
        iterations
      );

      if(iterations >= 1000) {
        throw std::logic_error("Could not find lambda to satisfy equation");
      }

      const FloatType lambda = (root_result.first + root_result.second) / 2;

      // Calculate uMin and calculate p from it
      const Vector2Type uMin = -(BTilde + lambda * BBar).inverse() * gTilde;
      return uMin(0) * g + uMin(1) * Binverse * g;
    }

    VectorType determineDirection(const FloatType trustRadius) const {
      if(positiveDefinite(hessians.current)) {
        return subspaceMinimize(trustRadius, gradients.current, hessians.current);
      }

      /* Now there are several possible cases for what the hessian is currently
       * like, but we'll need the eigenvalues to determine what to do.
       *
       * The hessian is a symmetric real matrix, so it is also self-adjoint,
       * which we exploit for faster eigenvalue calculation:
       */
      VectorType eigenvalues = hessians.current.template selfadjointView<Eigen::Lower>().eigenvalues();

      /* First case: there are negative eigenvalues. In this case we adjust
       * the hessian with a unitary matrix to become positive definite:
       */
      if((eigenvalues.array() < FloatType {0.0}).any()) {
        // Eigenvalues are ordered ascending
        const FloatType mostNegativeEigenvalue = eigenvalues(0);
        /* Nocedal & Wright say to choose alpha between (-λ,-2λ], so we
         * go right in the middle because we don't know anything about the
         * significance of the choice.
         */
        const FloatType alpha = FloatType {-1.5} * mostNegativeEigenvalue;
        const unsigned P = parameters.current.size();
        const Eigen::MatrixXd modifiedHessian = hessians.current + alpha * Eigen::MatrixXd::Identity(P, P);

        /* Subspace minimize the modified hessian
         *
         * NOTE: Nocedal & Wright say to differentiate some things further
         * here, but I think the subspace minimization procedure itself
         * differentiates those cases.
         */
        return subspaceMinimize(trustRadius, gradients.current, modifiedHessian);
      }

      /* Remaining case: There are no negative eigenvalues, but the hessian is
       * not positive definite. There must be eigenvalues of value zero. In
       * this case, we use the cauchy point to make progress:
       */
      const FloatType gradientNorm = gradients.current.norm();
      const FloatType tau = [&]() -> FloatType {
        const FloatType cauchyIntermediate = gradients.current.transpose() * hessians.current * gradients.current;
        if(cauchyIntermediate <= 0) {
          return 1.0;
        }

        return std::min(
          std::pow(gradientNorm, 3) / (trustRadius * cauchyIntermediate),
          1.0
        );
      }();

      return - (tau * trustRadius / gradientNorm) * gradients.current;
    }

    template<typename UpdateFunction>
    void prepare(UpdateFunction&& function, const FloatType trustRadius) {
      direction = determineDirection(trustRadius);
      parameters.proposed.noalias() = parameters.current + direction;
      function(parameters.proposed, values.proposed, gradients.proposed, hessians.proposed);
    }

    FloatType modelAgreement() const {
      const FloatType predictedValue = (
        values.current
        + detail::resolve(gradients.current.transpose() * direction)
        + detail::resolve(FloatType {0.5} * direction.transpose() * hessians.current * direction)
      );

      return (
        (values.current - values.proposed)
        / (values.current - predictedValue)
      );
    }

    template<typename UpdateFunction>
    void accept(UpdateFunction&& function, const FloatType trustRadius) {
      parameters.propagate();
      values.propagate();
      gradients.propagate();
      hessians.propagate();

      prepare(function, trustRadius);
    }
  };
};

} // namespace temple

#endif
