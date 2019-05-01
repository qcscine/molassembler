/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#ifndef INCLUDE_TEMPLE_LBFGS_H
#define INCLUDE_TEMPLE_LBFGS_H

#include <Eigen/Core>
#include <iostream>

namespace temple {

namespace detail {

struct DoNothingObserver {
  template<typename T>
  void operator() (T&& /* x */) {}
};

// This is similar to std::clamp except it doesn't expect lo <= hi (a <= b)
template<typename FloatType>
FloatType clamp(
  const FloatType a,
  const FloatType b,
  const FloatType value
) {
  if(a < b) {
    if(value < a) {
      return a;
    }

    if(b < value) {
      return b;
    }
  }

  if(value < b) {
    return b;
  }

  if(value > a) {
    return a;
  }

  return value;
}

template<typename FloatType>
FloatType polynomialMinimumExtrapolation(
  const FloatType f0,
  const FloatType d0,
  const FloatType f1,
  const FloatType d1,
  const FloatType limit = 1
) {
  const FloatType n = 3 * (f1 - f0) - 2 * d0 - d1;
  const FloatType e = d0 + d1 - 2 * (f1 - f0);

  // find the minimum of the derivative of the polynomial

  FloatType temp = std::max(n * n - 3 * e * d0, FloatType {0.0});

  if(temp < 0) {
    return 0.5;
  }

  temp = std::sqrt(temp);

  if(std::abs(e) <= std::numeric_limits<FloatType>::epsilon()) {
    return 0.5;
  }

  // figure out the two possible min values
  const FloatType x1 = (temp - n) / (3 * e);
  const FloatType x2 = -(temp + n) / (3 * e);

  // compute the value of the interpolating polynomial at these two points
  const FloatType y1 = f0 + d0 * x1 + n * x1 * x1 + e * x1 * x1 * x1;
  const FloatType y2 = f0 + d0 * x2 + n * x2 * x2 + e * x2 * x2 * x2;

  // pick the best point
  if(y1 < y2) {
    return clamp(FloatType {0}, limit, x1);
  }

  return clamp(FloatType {0}, limit, x2);
}

template<typename FloatType, typename UpdateFunction>
std::pair<FloatType, unsigned> lineSearch(
  UpdateFunction&& function,
  const Eigen::Ref<Eigen::Matrix<FloatType, Eigen::Dynamic, 1>>& parameters,
  const Eigen::Matrix<FloatType, Eigen::Dynamic, 1>& direction,
  const FloatType pivotValue,
  const FloatType pivotGradientDotDirection,
  const FloatType rho,
  const FloatType sigma,
  const FloatType goalValue,
  const unsigned maximumIterations
) {
  assert(0 < rho && rho < sigma && sigma < 1 && maximumIterations > 0);

  /* A note to the algorithm source:
   * The bracketing phase of this function is
   * implemented according to block 2.6.2 from the book Practical Methods of
   * Optimization by R. Fletcher.  The sectioning phase is an implementation of
   * 2.6.4 from the same book.
   */

  /* These parameters control the alpha jump size during the bracketing phase
   * of the search.
   */
  constexpr FloatType tau1a = 1.4;
  constexpr FloatType tau1b = 9;

  static_assert(
    FloatType {1} <= tau1a && tau1a < tau1b,
    "For this algorithm to function correctly, 1 <= tau1a < tau1b"
  );

  constexpr FloatType tau2 = 1.0 / 10.0;
  constexpr FloatType tau3 = 1.0 / 2.0;

  static_assert(
    FloatType {0} < tau2 && tau2 < tau3 && tau3 <= FloatType {0.5},
    "For this algorithm to function correctly, 0 < tau2 < tau3 <= 0.5 must be true."
  );

  /* This function should yield the number of iterations and an alpha value
   * but readability should not be detracted if possible.
   */
  unsigned long iteration = 0;
  auto constructResult = [&iteration](FloatType alphaValue) -> std::pair<FloatType, unsigned> {
    return std::make_pair(alphaValue, iteration);
  };

  /* Stop right away and return a step size of 0 if the gradient is 0 at the
   * starting point
   */
  if(std::abs(pivotGradientDotDirection) <= std::abs(pivotValue) * std::numeric_limits<FloatType>::epsilon()) {
    return constructResult(0);
  }

  // Stop right away if the current value is good enough according to goalValue
  if(pivotValue <= goalValue) {
    return constructResult(0);
  }

  // Figure out a reasonable upper bound on how large alpha can get
  const FloatType mu = (goalValue - pivotValue) / (rho * pivotGradientDotDirection);

  FloatType alpha = 1;
  if(mu < 0) {
    alpha = -alpha;
  }
  alpha = clamp(FloatType {0}, FloatType {0.65} * mu, alpha);

  FloatType lastAlpha = 0;
  FloatType lastValue = pivotValue;
  FloatType lastGradientDotDirectionValue = pivotGradientDotDirection;

  // The bracketing stage will find a range of points [a,b]
  // that contains a reasonable solution to the line search
  FloatType a, b;

  // These variables will hold the values and derivatives of f(a) and f(b)
  FloatType a_val, b_val, a_val_der, b_val_der;

  // This threshold value represents the Wolfe curvature condition
  const FloatType threshold = std::abs(sigma * pivotGradientDotDirection);

  // Gradient vector for storing results to objective function calls
  Eigen::VectorXd gradient(parameters.size());

  /* Bracketing stage: Find an interval [a,b] for alpha which is known to
   * contain acceptable continuation points
   */
  while (true) {
    ++iteration;
    FloatType value;
    function(parameters + direction * alpha, value, gradient);
    const FloatType gradientDotDirectionValue = gradient.dot(direction);

    // we are done with the line search since we found a value smaller
    // than the minimum f value
    if(value <= goalValue) {
      return constructResult(alpha);
    }

    if(
      value > pivotValue + rho * alpha * pivotGradientDotDirection
      || value >= lastValue
    ) {
      a_val = lastValue;
      a_val_der = lastGradientDotDirectionValue;
      b_val = value;
      b_val_der = gradientDotDirectionValue;

      a = lastAlpha;
      b = alpha;
      break;
    }

    if(std::abs(gradientDotDirectionValue) <= threshold) {
      return constructResult(alpha);
    }

    // if we are stuck not making progress then quit with the current alpha
    if(lastAlpha == alpha || iteration >= maximumIterations) {
      return constructResult(alpha);
    }

    if(gradientDotDirectionValue >= 0) {
      a_val = value;
      a_val_der = gradientDotDirectionValue;
      b_val = lastValue;
      b_val_der = lastGradientDotDirectionValue;

      a = alpha;
      b = lastAlpha;
      break;
    }

    const FloatType alphaBackup = alpha;

    // Pick a larger range [first, last] for alpha
    FloatType first, last;
    if(mu > 0) {
      first = std::min(mu, alpha + tau1a * (alpha - lastAlpha));
      last  = std::min(mu, alpha + tau1b * (alpha - lastAlpha));
    } else {
      first = std::max(mu, alpha + tau1a * (alpha - lastAlpha));
      last  = std::max(mu, alpha + tau1b * (alpha - lastAlpha));
    }

    // Pick a point between first and last by doing some kind of interpolation
    if(lastAlpha < alpha) {
      alpha = (
        lastAlpha
        + (alpha - lastAlpha) * polynomialMinimumExtrapolation<FloatType>(
          lastValue,
          lastGradientDotDirectionValue,
          value,
          gradientDotDirectionValue,
          1e10
        )
      );
    } else {
      alpha = (
        alpha
        + (lastAlpha - alpha) * polynomialMinimumExtrapolation<FloatType>(
          value,
          gradientDotDirectionValue,
          lastValue,
          lastGradientDotDirectionValue,
          1e10
        )
      );
    }

    alpha = clamp(first, last, alpha);

    lastAlpha = alphaBackup;

    lastValue = value;
    lastGradientDotDirectionValue = gradientDotDirectionValue;
  }

  /* Sectioning phase: Make the interval determined in the bracketing phase
   * smaller each iteration.
   *
   * Book reference: from section 2.6.4
   */
  while (true) {
    ++iteration;
    FloatType first = a + tau2 * (b - a);
    FloatType last = b - tau3 * (b - a);

    // use interpolation to pick alpha between first and last
    alpha = a + (b - a) * polynomialMinimumExtrapolation<FloatType>(a_val, a_val_der, b_val, b_val_der);
    alpha = clamp(first, last, alpha);

    FloatType value;
    function(parameters + alpha * direction, value, gradient);
    const FloatType gradientDotDirectionValue = gradient.dot(direction);

    /* We are done with the line search since we found a value smaller than the
     * minimum value or we ran out of iterations.
     */
    if(value <= goalValue || iteration >= maximumIterations) {
      return constructResult(alpha);
    }

    /* Stop if the interval gets so small that it isn't shrinking any more due
     * to rounding error
     */
    if(a == first || b == last) {
      return constructResult(b);
    }

    /* If alpha has basically become zero then just stop. Think of it like
     * this, if we take the largest possible alpha step will the objective
     * function change at all? If not then there isn't any point looking for a
     * better alpha.
     */
    const FloatType max_possible_alpha = std::max(std::abs(a),std::abs(b));
    if(
      std::abs(max_possible_alpha * pivotGradientDotDirection)
      <= std::abs(pivotValue) * std::numeric_limits<FloatType>::epsilon()
    ) {
      return constructResult(alpha);
    }

    if(
      value > pivotValue + rho * alpha * pivotGradientDotDirection
      || value >= a_val
    ) {
      b = alpha;
      b_val = value;
      b_val_der = gradientDotDirectionValue;
    } else {
      if(std::abs(gradientDotDirectionValue) <= threshold) {
        return constructResult(alpha);
      }

      if((b - a) * gradientDotDirectionValue >= 0) {
        b = a;
        b_val = a_val;
        b_val_der = a_val_der;
      }

      a = alpha;
      a_val = value;
      a_val_der = gradientDotDirectionValue;
    }
  }
}

} // namespace detail

template<typename FloatType, unsigned ringBufferSize>
struct LBFGS {
  constexpr static bool isPowerOfTwo(unsigned x) {
    return (x & (x - 1)) == 0;
  }

  static_assert(
    isPowerOfTwo(ringBufferSize) && ringBufferSize != 0,
    "For good performance, the ring buffer size has to be a power of two"
  );

  using MatrixType = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;

  template<typename UpdateFunction, typename Checker, typename Observer = detail::DoNothingObserver>
  unsigned optimize(
    Eigen::Ref<VectorType> parameters,
    UpdateFunction&& function,
    Checker&& check,
    Observer&& observer = Observer {}
  ) {
    /* number of parameters treated */
    unsigned nParams = parameters.size();
    FloatType value, oldValue;

    // Ring buffered old parameters and gradients
    unsigned ringBufferElementCount = 0, ringBufferOffset = 0;
    Eigen::Matrix<FloatType, Eigen::Dynamic, ringBufferSize> dg(nParams, ringBufferSize);
    Eigen::Matrix<FloatType, Eigen::Dynamic, ringBufferSize> dx(nParams, ringBufferSize);
    Eigen::Matrix<FloatType, ringBufferSize, 1> dxDotDg;

    // TODO these calls should not be necessary if the circular buffer work is correct
    dg.setZero();
    dx.setZero();
    dxDotDg.setZero();

    /* Setting all gradients to zero. */
    VectorType gradients = VectorType::Zero(nParams);
    /* Setting all old gradients to zero. */
    VectorType oldGradients = VectorType::Zero(nParams);
    /* Set "assumed previous parameters" as oldParameters. */
    VectorType oldParameters(parameters);

    // Calculate initial function value and gradients
    function(parameters, value, gradients);

    /* start with one steepest descent step */
    VectorType stepVector = FloatType {-0.1} * stepLength * gradients;

    unsigned cycle = 0;
    while (true) {
      ++cycle;

      // Store old results in local variables for possible backtracking needs
      oldParameters = parameters;
      oldGradients = gradients;
      oldValue = value;

      /* Update the parameters by the direction multiplied by stepLength */
      parameters.noalias() += (stepLength * stepVector);
      /* Update gradients/value and check convergence */
      function(parameters, value, gradients);

      observer(parameters);

      if(check.checkMaxIterations(cycle) || check.checkConvergence(parameters, value, gradients)) {
        break;
      }

      /* Generate a new direction */

      // Armijo condition (Wolfe condition I), used keep the step length sufficient
      bool armijo = value <= (oldValue + c1 * stepLength * oldGradients.dot(stepVector));
      /* Curvature condition (Wolfe condition II) */
      bool curvature = gradients.dot(stepVector) >= (c2 * oldGradients.dot(stepVector));
      /* Backtracking condition */
      bool backtracking = value > oldValue;

      /* Linesearch */
      if(linesearch) {
        if(!armijo) {
          /* Check need for backtracking */
          if(backtracking) {
            //std::cout << "Backtracking!\n";
            parameters = oldParameters;
            gradients = oldGradients;
            /* Adjust step length */
            stepLength *= FloatType {0.5};
            continue;
          }

          /* Technically the step size is wrong and there should be backtracking
           * followed by a step with the correct length either way.
           * However, if the step was sane (i.e. the value fell), the step
           * is kept in order to save computation time on updates.
           */
          /* Adjust step length */
          stepLength *= FloatType {0.5};
        } else if(!curvature) {
          /* Check need for backtracking */
          if(backtracking) {
            //std::cout << "Backtracking!\n";
            parameters = oldParameters;
            gradients = oldGradients;
            /* Adjust step length */
            stepLength *= FloatType {1.5};
            continue;
          }

          /* Technically the step size is wrong and there should be backtracking
           * followed by a step with the correct length either way.
           * However, if the step was sane (i.e. the value fell), the step
           * is kept in order to save computation time on updates.
           */
          /* Adjust step length */
          stepLength *= FloatType {1.5};
        }
      }

      /* Add current gradients and parameters to storage */
      if(ringBufferElementCount < ringBufferSize) {
        dg.col(ringBufferElementCount) = gradients - oldGradients;
        dx.col(ringBufferElementCount) = parameters - oldParameters;
        dxDotDg(ringBufferElementCount) = dx.col(ringBufferElementCount).dot(dg.col(ringBufferElementCount));
        ++ringBufferElementCount;
      } else {
        // Rotate columns without copying (introduce modulo offset)
        const unsigned columnOffset = (ringBufferElementCount + ringBufferOffset) % ringBufferSize;
        dg.col(columnOffset) = gradients - oldGradients;
        dx.col(columnOffset) = parameters - oldParameters;
        dxDotDg(columnOffset) = dx.col(columnOffset).dot(dg.col(columnOffset));
        ringBufferOffset = (ringBufferOffset + 1) % ringBufferSize;
      }

      /* L-BFGS update */
      VectorType alpha(ringBufferElementCount);
      stepVector = -gradients;
      auto newestToOldest = [&](unsigned i) {
        FloatType dxDotdg = dxDotDg(i);
        if(std::fabs(dxDotdg) < FloatType {1.0e-6}) {
          alpha[i] = dx.col(i).dot(stepVector) / FloatType {1.0e-6};
        } else {
          alpha[i] = dx.col(i).dot(stepVector) / dxDotdg;
        }

        stepVector.noalias() -= alpha[i] * dg.col(i);
      };

      // Newest to oldest
      const int newestOffset = (ringBufferElementCount + ringBufferOffset - 1) % ringBufferSize;
      for(int i = newestOffset; i > - 1; --i) {
        newestToOldest(i);
      }
      for(
        int i = std::min(ringBufferElementCount - 1, ringBufferSize - 1);
        i > newestOffset;
        --i
      ) {
        newestToOldest(i);
      }

      stepVector *= dxDotDg(newestOffset) / dg.col(newestOffset).squaredNorm();

      auto oldestToNewest = [&](unsigned i) {
        FloatType dxDotdg = dxDotDg(i);
        FloatType beta = dg.col(i).dot(stepVector);

        if(std::fabs(dxDotdg) < FloatType {1.0e-6}) {
          beta /= FloatType {1.0e-6};
        } else {
          beta /= dxDotDg(i);
        }

        stepVector.noalias() += (alpha[i] - beta) * dx.col(i);
      };

      // Oldest to newest
      const unsigned oldestOffset = (ringBufferElementCount + ringBufferOffset) % ringBufferSize;
      for(unsigned i = oldestOffset; i < ringBufferElementCount; ++i) {
        oldestToNewest(i);
      }
      for(unsigned i = 0; i < oldestOffset; ++i) {
        oldestToNewest(i);
      }

      /*unsigned additionalCycles;
      std::tie(stepLength, additionalCycles) = detail::lineSearch(
        function,
        parameters,
        stepVector,
        value,
        gradients.dot(stepVector),
        c1,
        c2,
        FloatType {0.0},
        100u
      );
      cycle += additionalCycles;*/
      //std::cout << stepLength << " " << additionalCycles << "\n";
      //std::cout << stepLength << "\n";
    }

    return cycle;
  }

  /// @brief Switch to turn on and off the use of a linesearch
  bool linesearch = true;
  /**
   * @brief 1st parameter for the Wolfe conditions.
   *
   * This parameter is only relevant if the linesearch is turned on.
   */
  FloatType c1 = 0.0001;
  /**
   * @brief 2nd parameter for the Wolfe conditions.
   *
   * This parameter is only relevant if the linesearch is turned on.
   * Also a value as low as 0.1 can be used.
   */
  FloatType c2 = 0.9;
  /**
   * @brief The initial step length used in the L-BFGS.
   *
   * Note: the first step is a gradient descent with 0.1 times the steplength.
   */
  FloatType stepLength = 1.0;
};

} // namespace temple

#endif
