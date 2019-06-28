/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#ifndef INCLUDE_TEMPLE_LBFGS_H
#define INCLUDE_TEMPLE_LBFGS_H

#include <Eigen/Core>
#include <iostream>
#include "temple/STL17.h"

/* TODO
 * - Detect ill-conditioning of LBFGS? Doesn't do well at all near maxima
 */

namespace temple {

namespace detail {

template<typename FloatType> constexpr FloatType cfabs(FloatType x) noexcept {
  return (x >= 0) ? x : -x;
}

struct DoNothingObserver {
  template<typename T>
  void operator() (T&& /* x */) {}
};

template<typename FloatType>
struct FloatUpdateBuffer {
  static_assert(
    std::is_floating_point<FloatType>::value,
    "This struct is not intended for anything else but floats"
  );

  FloatType old;
  FloatType current;

  //! Returns signed difference: current - old
  FloatType absDelta() const {
    return std::fabs(current - old);
  }

  void propagate() { old = current; }
};

template<typename EigenType>
struct EigenUpdateBuffer {
  EigenType old;
  EigenType current;

  typename EigenType::Scalar deltaNorm() const {
    return (current - old).norm();
  }

  void propagate() {
    /* Swapping instead of moving current into old is three moves instead of
     * one, but avoids allocation of new memory (after moving from current, we
     * need to resize current to have the same size). Assuming moves are
     * cheap (trade ownership of allocated memory) and allocation could be
     * costly, this might be cheaper.
     */
    std::swap(old, current);
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

template<typename FloatType, typename VectorType, typename UpdateFunction>
struct ClampFunction {
  UpdateFunction function;
  const VectorType& low;
  const VectorType& high;

  ClampFunction(
    UpdateFunction&& fn,
    const VectorType& lo,
    const VectorType& hi
  ) : function(fn), low(lo), high(hi) {}

  void operator() (
    const VectorType& parameters,
    FloatType& value,
    Eigen::Ref<VectorType> gradient
  ) {
    function(
      parameters.cwiseMax(low).cwiseMin(high),
      value,
      gradient
    );
  }
};

template<typename VectorType, typename UpdateFunction>
auto clampFunction(
  UpdateFunction&& fn,
  const VectorType& lo,
  const VectorType& hi
) {
  using FloatType = typename VectorType::Scalar;
  return ClampFunction<FloatType, VectorType, UpdateFunction>(
    std::forward<UpdateFunction>(fn),
    lo,
    hi
  );
}

namespace boxes {

template<typename VectorType>
void clampToBounds(Eigen::Ref<VectorType> /* parameters */) {}

template<typename VectorType, typename BoxType>
void clampToBounds(Eigen::Ref<VectorType> parameters, const BoxType& box) {
  box.clamp(parameters);
}


template<typename VectorType>
bool validate(const VectorType& /* parameters */) {
  return true;
}

template<typename VectorType, typename BoxType>
bool validate(const VectorType& parameters, const BoxType& box) {
  return box.validate(parameters);
}

template<typename VectorType>
void adjustGradient(
  Eigen::Ref<VectorType> /* gradient */,
  const VectorType& /* parameters */
) {}

template<typename VectorType, typename BoxType>
void adjustGradient(
  Eigen::Ref<VectorType> gradient,
  const VectorType& parameters,
  const BoxType& box
) {
  box.adjustGradient(gradient, parameters);
}

template<typename VectorType>
void adjustDirection(
  Eigen::Ref<VectorType> /* direction */,
  const VectorType& /* parameters */,
  const VectorType& /* gradient */
) {}

template<typename VectorType, typename BoxType>
void adjustDirection(
  Eigen::Ref<VectorType> direction,
  const VectorType& parameters,
  const VectorType& gradient,
  const BoxType& box
) {
  box.adjustDirection(direction, parameters, gradient);
}

} // namespace boxes

} // namespace detail

template<typename FloatType = double, unsigned ringBufferSize = 16>
class LBFGS {
public:
  constexpr static bool isPowerOfTwo(unsigned x) {
    return (x & (x - 1)) == 0;
  }

  static_assert(
    isPowerOfTwo(ringBufferSize) && ringBufferSize != 0,
    "For good performance, the ring buffer size has to be a power of two"
  );

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

  /**
   * @brief Type storing boundaries for parameter values
   */
  struct Box {
    //! Lower bounds for parameters
    VectorType minima;
    //! Upper bounds for parameters
    VectorType maxima;

    Box() = default;
    Box(const VectorType& passMinima, const VectorType& passMaxima)
      : minima(passMinima),
        maxima(passMaxima)
    {
      assert((minima.array() <= maxima.array()).all());
    }

    //! Clamp parameter values to the bounds
    void clamp(Eigen::Ref<VectorType> x) const {
      x = x.cwiseMax(minima).cwiseMin(maxima);
    }

    //! Test whether parameters are within the bounds
    bool validate(const VectorType& v) const {
      return (
        (minima.array() <= v.array()).all()
        && (v.array() <= maxima.array()).all()
      );
    }

    //! Zero out those elements of a gradient which are bounded
    void adjustGradient(
      Eigen::Ref<VectorType> gradient,
      const VectorType& parameters
    ) const {
      constexpr FloatType gapEpsilon = 1e-8;
      const Eigen::Index P = parameters.size();

      for(Eigen::Index i = 0; i < P; ++i) {
        if(
          (
            minima[i] + gapEpsilon >= parameters[i]
            && gradient[i] > 0
          ) || (
            maxima[i] - gapEpsilon <= parameters[i]
            && gradient[i] < 0
          )
        ) {
          gradient[i] = 0;
        }
      }
    }

    //! Adjust a direction vector to avoid exceeding the box
    void adjustDirection(
      Eigen::Ref<VectorType> direction,
      const VectorType& parameters,
      const VectorType& gradient
    ) const {
      constexpr FloatType gapEpsilon = 1e-8;

      const Eigen::Index P = parameters.size();
      for(Eigen::Index i = 0; i < P; ++i) {
        if(
          minima[i] + gapEpsilon >= parameters[i]
          && gradient[i] > 0
        ) {
          direction[i] = minima[i] - parameters[i];
        } else if(
          maxima[i] - gapEpsilon <= parameters[i]
          && gradient[i] < 0
        ) {
          direction[i] = maxima[i] - parameters[i];
        }
      }
    }
  };

  /**
   * @brief Class carrying proposed parameter updates in an optimization
   */
  struct StepValues {
    detail::EigenUpdateBuffer<VectorType> parameters;
    detail::FloatUpdateBuffer<FloatType> values;
    detail::EigenUpdateBuffer<VectorType> gradients;

    /**
     * @brief Initialize class state and generate first direction vector
     *
     * @param function objective function generating values and gradients from parameters
     * @param initialParameters Parameters passed to optimization method
     * @param multiplier Initial step length
     *
     * @return First direction vector
     */
    template<typename UpdateFunction, typename ... Boxes>
    VectorType generateInitialDirection(
      UpdateFunction&& function,
      Eigen::Ref<VectorType> initialParameters,
      const FloatType multiplier,
      const Boxes& ... boxes
    ) {
      const unsigned P = initialParameters.size();
      parameters.old = std::move(initialParameters);
      gradients.old.resize(P);
      gradients.current.resize(P);

      // Initialize old values and gradients
      function(parameters.old, values.old, gradients.old);
      detail::boxes::adjustGradient<VectorType>(gradients.old, parameters.old, boxes ...);

      // Initialize new with a small steepest descent step
      VectorType direction;
      const FloatType gradientNorm = gradients.old.norm();
      if(gradientNorm > FloatType {1e-4}) {
        direction = -gradients.old / gradientNorm;
      } else {
        direction = FloatType {-1.0} * gradients.old;
      }

      prepare(function, multiplier, direction, boxes ...);

      return direction;
    }

    template<typename UpdateFunction, typename ... Boxes>
    void prepare(
      UpdateFunction&& function,
      const FloatType multiplier,
      const VectorType& direction,
      const Boxes& ... boxes
    ) {
      parameters.current.noalias() = parameters.old + multiplier * direction;
      detail::boxes::clampToBounds<VectorType>(parameters.current, boxes ...);
      function(parameters.current, values.current, gradients.current);
      detail::boxes::adjustGradient<VectorType>(gradients.current, parameters.current, boxes ...);
    }

    /**
     * @brief Accept a step, create new parameters and evaluate the function
     *   to get value and gradients
     *
     * @param function objective function generating values and gradients from parameters
     * @param multiplier Step length
     * @param direction Direction to advance parameters
     */
    template<typename UpdateFunction, typename ... Boxes>
    void propagate(
      UpdateFunction&& function,
      const FloatType multiplier,
      const VectorType& direction,
      const Boxes& ... boxes
    ) {
      parameters.propagate();
      values.propagate();
      gradients.propagate();

      prepare(function, multiplier, direction, boxes ...);
    }
  };


  /**
   * @brief Find parameters to minimize an objective function
   *
   * @tparam UpdateFunction A callable object taking arguments (const VectorType&,
   *   double&, Ref<VectorType>). Any return values are ignored.
   * @tparam Checker An object implementing two methods:
   *   - checkMaxIterations: unsigned -> bool
   *   - checkConvergence: (const VectorType&, double, const VectorType&) -> bool
   *     The arguments passed are, in order: parameters, value and gradient
   *   These functions should return true if optimization should cease
   * @tparam Observer An observer callable object taking parameters. Return
   *   values are ignored.
   * @param parameters The initial parameters to optimize. Resulting
   *   optimization parameters are written to this argument.
   * @param function The objective function to optimize.
   * @param check The Checker instance
   * @param observer The Observer instance
   *
   * @return A @p OptimizationReturnType instance
   */
  template<
    typename UpdateFunction,
    typename Checker,
    typename Observer = detail::DoNothingObserver
  > OptimizationReturnType minimize(
    Eigen::Ref<VectorType> parameters,
    UpdateFunction&& function,
    Checker&& check,
    Observer&& observer = Observer {}
  ) {
    return minimizeBase(
      parameters,
      std::forward<UpdateFunction>(function),
      std::forward<Checker>(check),
      std::forward<Observer>(observer)
    );
  }

  /**
   * @brief Find parameters to maximize an objective function
   *
   * @tparam UpdateFunction A callable object taking arguments (const VectorType&,
   *   double&, Ref<VectorType>). Any return values are ignored.
   * @tparam Checker An object implementing two methods:
   *   - checkMaxIterations: unsigned -> bool
   *   - checkConvergence: (const VectorType&, double, const VectorType&) -> bool
   *     The arguments passed are, in order: parameters, value and gradient
   *   These functions should return true if optimization should cease
   * @tparam Observer An observer callable object taking parameters. Return
   *   values are ignored.
   * @param parameters The initial parameters to optimize. Resulting
   *   optimization parameters are written to this argument.
   * @param function The objective function to optimize.
   * @param check The Checker instance
   * @param observer The Observer instance
   *
   * @return The number of calls to @p function made during optimization
   */
  template<
    typename UpdateFunction,
    typename Checker,
    typename Observer = detail::DoNothingObserver
  > OptimizationReturnType maximize(
    Eigen::Ref<VectorType> parameters,
    UpdateFunction&& function,
    Checker&& check,
    Observer&& observer = Observer {}
  ) {
    OptimizationReturnType results = minimize(
      parameters,
      detail::negateFunction<VectorType>(
        std::forward<UpdateFunction>(function)
      ),
      std::forward<Checker>(check),
      std::forward<Observer>(observer)
    );

    // Revert sign of negated function minimization
    results.value *= -1;
    results.gradient *= -1;

    return results;
  }

  /**
   * @brief Find parameters to minimize an objective function within bounds
   *
   * @tparam UpdateFunction A callable object taking arguments (const VectorType&,
   *   double&, Ref<VectorType>). Any return values are ignored.
   * @tparam Checker An object implementing two methods:
   *   - checkMaxIterations: unsigned -> bool
   *   - checkConvergence: (const VectorType&, double, const VectorType&) -> bool
   *     The arguments passed are, in order: parameters, value and gradient
   *   These functions should return true if optimization should cease
   * @tparam Observer An observer callable object taking parameters. Return
   *   values are ignored.
   * @param parameters The initial parameters to optimize. Resulting
   *   optimization parameters are written to this argument.
   * @param box Bounds on the parameters
   * @param function The objective function to optimize.
   * @param check The Checker instance
   * @param observer The Observer instance
   *
   * @note Gradients passed to @p check do not contain bounded components
   *
   * @return The number of calls to @p function made during optimization
   */
  template<
    typename UpdateFunction,
    typename Checker,
    typename Observer = detail::DoNothingObserver
  > OptimizationReturnType minimize(
    Eigen::Ref<VectorType> parameters,
    const Box& box,
    UpdateFunction&& function,
    Checker&& check,
    Observer&& observer = Observer {}
  ) {
    return minimizeBase(
      parameters,
      std::forward<UpdateFunction>(function),
      std::forward<Checker>(check),
      std::forward<Observer>(observer),
      box
    );
  }

  /**
   * @brief Find parameters to maximize an objective function within bounds
   *
   * @tparam UpdateFunction A callable object taking arguments (const VectorType&,
   *   double&, Ref<VectorType>). Any return values are ignored.
   * @tparam Checker An object implementing two methods:
   *   - checkMaxIterations: unsigned -> bool
   *   - checkConvergence: (const VectorType&, double, const VectorType&) -> bool
   *     The arguments passed are, in order: parameters, value and gradient
   *   These functions should return true if optimization should cease
   * @tparam Observer An observer callable object taking parameters. Return
   *   values are ignored.
   * @param parameters The initial parameters to optimize. Resulting
   *   optimization parameters are written to this argument.
   * @param box Bounds on the parameters
   * @param function The objective function to optimize.
   * @param check The Checker instance
   * @param observer The Observer instance
   *
   * @return The number of calls to @p function made during optimization
   */
  template<
    typename UpdateFunction,
    typename Checker,
    typename Observer = detail::DoNothingObserver
  > OptimizationReturnType maximize(
    Eigen::Ref<VectorType> parameters,
    const Box& box,
    UpdateFunction&& function,
    Checker&& check,
    Observer&& observer = Observer {}
  ) {
    OptimizationReturnType results = minimize(
      parameters,
      box,
      detail::negateFunction<VectorType>(
        std::forward<UpdateFunction>(function)
      ),
      std::forward<Checker>(check),
      std::forward<Observer>(observer)
    );

    // Revert sign of negated function minimization
    results.value *= -1;
    results.gradient *= -1;

    return results;
  }

  /**
   * @brief 1st parameter for the Wolfe conditions.
   */
  FloatType c1 = 0.0001;
  /**
   * @brief 2nd parameter for the Wolfe conditions.
   *
   * A value as low as 0.1 can be used.
   */
  FloatType c2 = 0.9;
  /**
   * @brief The initial step length used in the L-BFGS.
   *
   * Note: the first step is a gradient descent with 0.1 times the steplength.
   */
  FloatType stepLength = 1.0;

private:
  /**
   * @brief Class containing hessian information for LBFGS update
   */
  struct CollectiveRingBuffer {
    //! Ring buffer containing gradient differences: y_k = g_{k+1} - g_k
    Eigen::Matrix<FloatType, Eigen::Dynamic, ringBufferSize> y;
    //! Ring buffer containing parameter differnces: s_k = x_{k+1} - x_k
    Eigen::Matrix<FloatType, Eigen::Dynamic, ringBufferSize> s;
    //! Ring buffered result of s_k.dot(y_k)
    Eigen::Matrix<FloatType, ringBufferSize, 1> sDotY;
    unsigned count = 0, offset = 0;

    CollectiveRingBuffer(const unsigned nParams)
      : y(nParams, ringBufferSize),
        s(nParams, ringBufferSize)
    {}

    unsigned newestOffset() const {
      return (count + offset - 1) % ringBufferSize;
    }

    unsigned oldestOffset() const {
      return (count + offset) % ringBufferSize;
    }

    /**
     * @brief Calls a unary function with the indices of the ring buffer in
     *   sequence from newest to oldest entries
     *
     * @tparam UnaryFn A unary callable with argument convertible from int
     * @param fn The function to call
     */
    template<typename UnaryFn>
    void newestToOldest(UnaryFn&& fn) const {
      assert(count > 0);
      const int newest = newestOffset();
      for(int i = newest; i > - 1; --i) {
        fn(i);
      }

      const int countOrSizeLimit = std::min(count - 1, ringBufferSize - 1);
      for(int i = countOrSizeLimit; i > newest; --i) {
        fn(i);
      }
    }

    /**
     * @brief Calls a unary function with the indices of the ring buffer in
     *   sequence from oldest to newest entries
     *
     * @tparam UnaryFn A unary callable with argument convertible from int
     * @param fn The function to call
     */
    template<typename UnaryFn>
    void oldestToNewest(UnaryFn&& fn) const {
      assert(count > 0);
      const unsigned oldest = oldestOffset();
      for(unsigned i = oldest; i < count; ++i) {
        fn(i);
      }

      for(unsigned i = 0; i < oldest; ++i) {
        fn(i);
      }
    }

    bool addInformation(
      const detail::EigenUpdateBuffer<VectorType>& parameterBuffer,
      const detail::EigenUpdateBuffer<VectorType>& gradientBuffer
    ) {
      bool dotProductNotZero;
      if(count < ringBufferSize) {
        y.col(count).noalias() = gradientBuffer.current - gradientBuffer.old;
        s.col(count).noalias() = parameterBuffer.current - parameterBuffer.old;
        sDotY(count) = s.col(count).dot(y.col(count));
        dotProductNotZero = (sDotY(count) != 0);
        ++count;
      } else {
        // Rotate columns without copying (introduce modulo offset)
        const unsigned columnOffset = (count + offset) % ringBufferSize;
        y.col(columnOffset).noalias() = gradientBuffer.current - gradientBuffer.old;
        s.col(columnOffset).noalias() = parameterBuffer.current - parameterBuffer.old;
        sDotY(columnOffset) = s.col(columnOffset).dot(y.col(columnOffset));
        dotProductNotZero = (sDotY(columnOffset) != 0);
        offset = (offset + 1) % ringBufferSize;
      }

      return dotProductNotZero;
    }

    /**
     * @brief Use stored gradient and parameter differences to generate a new
     *   direction
     *
     * @param q The vector in which to store the direction
     * @param currentGradient The current gradient
     */
    void generateNewDirection(
      Eigen::Ref<VectorType> q,
      const VectorType& currentGradient
    ) const {
      /* Note: Naming follows Numerical Optimization (2006)'s naming scheme,
       * p.178.
       */
      const unsigned kMinusOne = newestOffset();

      q = -currentGradient;
      VectorType alpha(count);

      newestToOldest(
        [&](const unsigned i) {
          alpha[i] = s.col(i).dot(q) / sDotY(i);
          q.noalias() -= alpha[i] * y.col(i);
        }
      );

      /* gamma_k = s_{k-1}^T y_{k-1} / (y_{k-1}^T y_{k-1})
       *         = s_{k-1}.dot(y_{k-1}) / y_{k-1}.squaredNorm()
       *
       * H_k^0 = y_k * I (where I is the identity matrix for n dimensions)
       * r = H_k^0 * q (where q is the current direction)
       *
       * q is not needed again later in the algorithm, so instead of defining
       * a new vector r we just reuse q.
       *
       * -> q *= gamma_k;
       */
      q *= sDotY(kMinusOne) / y.col(kMinusOne).squaredNorm();

      oldestToNewest(
        [&](const unsigned i) {
          const FloatType beta = y.col(i).dot(q) / sDotY(i);
          q.noalias() += (alpha[i] - beta) * s.col(i);
        }
      );
    }

    template<typename ... Boxes>
    bool updateAndGenerateNewDirection(
      Eigen::Ref<VectorType> direction,
      const StepValues& step,
      const Boxes& ... boxes
    ) {
      /* L-BFGS update: Use old parameters and gradients to improve steepest
       * descent direction
       */
      if(!addInformation(step.parameters, step.gradients)) {
        return false;
      }

      generateNewDirection(direction, step.gradients.current);
      detail::boxes::adjustDirection(
        direction,
        step.parameters.current,
        step.gradients.current,
        boxes ...
      );

      return true;
    }
  };

  /**
   * @brief Adjusts @p stepLength member and ensures the current parameter
   *   adjustment lowers the function value.
   *
   * The chosen direction and step length may overshoot or undershoot, and the
   * function value might increase. We try to adapt the step length in order to
   * ensure the function value decreases and step length is more or less
   * optimal.
   *
   * @param function Objective function generating value and gradients from parameters
   * @param observer Observer called with parameter updates
   * @param step Suggested step in parameters, values and gradients
   * @param direction Suggested step direction
   *
   * @return Number of function evaluations made in adjusting the step
   */
  template<
    typename UpdateFunction,
    typename Observer,
    typename ... Boxes
  > unsigned adjustStepAlongDirection(
    UpdateFunction&& function,
    Observer&& observer,
    StepValues& step,
    const VectorType& direction,
    const Boxes& ... boxes
  ) {
    unsigned iterations = 0;

    /* Numerical optimization p. 178, paraphrased: Generated search direction
     * is typically well scaled, so that step length = 1 should be accepted in
     * most iterations.
     *
     * Therefore, we only allow buildup of lengthening of the step length, but
     * not shortening.
     */
    if(stepLength < 1) {
      stepLength = 1;
    }

    for(iterations = 0; stepLength > 0; ++iterations) {
      const FloatType oldGradDotDirection = step.gradients.old.dot(direction);

      /* Armijo rule (Wolfe condition I) */
      bool armijoRule = (
        step.values.current
        <= step.values.old + c1 * stepLength * oldGradDotDirection
      );

      /* Curvature condition (Wolfe condition II) */
      bool curvatureCondition = (
        step.gradients.current.dot(direction)
        >= c2 * oldGradDotDirection
      );

      /* Backtracking condition */
      bool backtrack = (step.values.current > step.values.old);

      // Decide whether to shorten or extend the step length
      constexpr FloatType shortenFactor = 0.5;
      constexpr FloatType lengthenFactor = 1.5;
      static_assert(
        detail::cfabs(shortenFactor * lengthenFactor - 1.0) > 0.1,
        "Shortening and lengthening factor should not multiply to give "
        "approximately one, this can cause oscillation"
      );

      if(!armijoRule) {
        stepLength *= shortenFactor;
      } else if(!curvatureCondition) {
        stepLength *= lengthenFactor;
      }

      /* If we don't have to backtrack, then don't re-evaluate along this
       * direction, but accept the step and perform the next evaluation at
       * the next position. We have improved the function value with the step
       * and possibly adjusted the step length using the Wolfe conditions,
       * that's enough.
       */
      if(backtrack && (!armijoRule || !curvatureCondition)) {
        // Retry shortened/lengthened step (no propagation needed)
        step.prepare(function, stepLength, direction, boxes ...);
        observer(step.parameters.current);

        // Handle encountered parameter boundaries
        if(
          sizeof...(boxes) > 0
          && (step.parameters.current - step.parameters.old).squaredNorm() < FloatType {1e-8}
        ) {
          break;
        }
      } else {
        break;
      }
    }

    if(stepLength == 0) {
      throw std::logic_error("Could not find a step length that reduces objective function value.");
    }

    return iterations;
  }

  template<
    typename UpdateFunction,
    typename Checker,
    typename Observer,
    typename ... Boxes
  > OptimizationReturnType minimizeBase(
    Eigen::Ref<VectorType> parameters,
    UpdateFunction&& function,
    Checker&& check,
    Observer&& observer,
    Boxes&& ... boxes
  ) {
    // If there is a box, make sure the parameters are valid
    assert(detail::boxes::validate(parameters, boxes ...));

    // Set up a first small conjugate gradient step
    StepValues step;
    VectorType direction = step.generateInitialDirection(function, parameters, stepLength, boxes ...);
    observer(step.parameters.current);
    //std::cout << "Propose step to " << step.parameters.current.transpose() << ", value " << step.values.current << "\n";

    /* Set up ring buffer to keep changes in gradient and parameters to
     * approximate the inverse Hessian with
     */
    CollectiveRingBuffer ringBuffer {static_cast<unsigned>(parameters.size())};

    // Begin optimization loop
    unsigned iteration;
    for(
      iteration = 1;
      check.shouldContinue(iteration, stl17::as_const(step));
      ++iteration
    ) {
      // Line search the chosen direction
      iteration += adjustStepAlongDirection(function, observer, step, direction, boxes ...);

      /* Add more information to approximate the inverse Hessian and use it to
       * generate a new direction.
       */
      if(!ringBuffer.updateAndGenerateNewDirection(direction, step, boxes ...)) {
        break;
      }

      /* Accept the current step and prepare a new prospective step using the
       * updated direction vector
       */
      //std::cout << "Accepting step to " << step.parameters.current.transpose() << "\n";
      step.propagate(function, stepLength, direction, boxes ...);
      //std::cout << "Propose step to " << step.parameters.current.transpose() << ", value " << step.values.current << "\n";
      observer(step.parameters.current);
    }

    // Copy the optimal parameters back into the in/out argument
    parameters = std::move(step.parameters.current);
    return {
      iteration,
      step.values.current,
      std::move(step.gradients.current)
    };
  }
};

} // namespace temple

#endif
