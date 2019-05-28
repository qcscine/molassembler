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

template<typename FloatType>
struct FloatUpdateBuffer {
  static_assert(
    std::is_floating_point<FloatType>::value,
    "This struct is not intended for anything else but floats"
  );

  FloatType old;
  FloatType current;

  void propagate() {
    old = current;
  }

  void revert() {
    current = old;
  }
};

template<typename EigenType>
struct EigenUpdateBuffer {
  EigenType old;
  EigenType current;

  void propagate() {
    /* Swapping instead of moving current into old is three moves instead of
     * one, but avoids allocation of new memory (after moving from current, we
     * need to resize current to have the same size). Assuming moves are
     * cheap (trade ownership of allocated memory) and allocation with context
     * switches could be costly, this might be cheaper.
     */
    std::swap(old, current);
  }

  void revert() {
    current = old;
  }
};

} // namespace detail

template<typename FloatType, unsigned ringBufferSize>
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

  struct StepValues {
    detail::EigenUpdateBuffer<VectorType> parameters;
    detail::EigenUpdateBuffer<FloatType> values;
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
    template<typename UpdateFunction>
    VectorType generateInitialDirection(
      UpdateFunction&& function,
      Eigen::Ref<VectorType> initialParameters,
      FloatType multiplier
    ) {
      const unsigned P = initialParameters.size();
      parameters.old = std::move(initialParameters);
      gradients.old.resize(P);
      gradients.current.resize(P);

      // Initialize old values and gradients
      function(parameters.old, values.old, gradients.old);

      // Initialize new with a small steepest descent step
      VectorType direction = FloatType {-0.1} * gradients.old;
      // Create new parameters
      parameters.current.noalias() = parameters.old + multiplier * direction;
      // Re-evaluate the function, assigning current value and gradient
      function(parameters.current, values.current, gradients.current);

      return direction;
    }

    template<typename UpdateFunction>
    void propagate(
      UpdateFunction&& function,
      const FloatType multiplier,
      const VectorType& direction
    ) {
      parameters.propagate();
      values.propagate();
      gradients.propagate();

      parameters.current.noalias() = parameters.old + multiplier * direction;
      function(parameters.current, values.current, gradients.current);
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
    typename Observer
  > unsigned adjustStepAlongDirection(
    UpdateFunction&& function,
    Observer&& observer,
    StepValues& step,
    Eigen::Ref<VectorType> direction
  ) {
    unsigned iterations = 0;

    for(iterations = 0; true; ++iterations) {
      /* Armijo condition (Wolfe condition I) */
      bool armijoCondition = (
        step.values.current
        <= step.values.old + c1 * stepLength * step.gradients.old.dot(direction)
      );

      /* Curvature condition (Wolfe condition II) */
      bool curvatureCondition = (
        step.gradients.current.dot(direction)
        >= c2 * step.gradients.old.dot(direction)
      );

      /* Backtracking condition */
      bool backtrack = (step.values.current > step.values.old);

      // Decide whether to shorten or extend the step length
      constexpr FloatType shortenFactor = 0.5;
      constexpr FloatType lengthenFactor = 1.5;
      static_assert(
        std::fabs(shortenFactor * lengthenFactor - 1.0) > 0.1,
        "Shortening and lengthening factor should not multiply to give "
        "approximately one, this can cause oscillation"
      );

      if(!armijoCondition) {
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
      if(backtrack && (!armijoCondition || !curvatureCondition)) {
        // Retry shortened/lengthened step (no propagation needed)
        step.parameters.current.noalias() = step.parameters.old + stepLength * direction;
        /* Update gradients and value */
        function(step.parameters.current, step.values.current, step.gradients.current);
        observer(step.parameters.current);
      } else {
        break;
      }
    }

    return iterations;
  }

  template<
    typename UpdateFunction,
    typename Checker,
    typename Observer = detail::DoNothingObserver
  > unsigned optimize(
    Eigen::Ref<VectorType> initialParameters,
    UpdateFunction&& function,
    Checker&& check,
    Observer&& observer = Observer {}
  ) {
    /* number of initialParameters treated */
    // Note: Do not reorder, StepValues moves from initialParameters
    const unsigned P = initialParameters.size();
    StepValues step;
    VectorType direction = step.generateInitialDirection(function, initialParameters, stepLength);
    CollectiveRingBuffer ringBuffer {P};
    observer(step.parameters.current);

    unsigned iteration;
    for(
      iteration = 1;
      (
        !check.checkMaxIterations(iteration)
        && !check.checkConvergence(step.parameters.current, step.values.current, step.gradients.current)
      );
      ++iteration
    ) {
      iteration += adjustStepAlongDirection(function, observer, step, direction);

      ringBuffer.updateAndGenerateNewDirection(direction, step);

      /* Accept the current step and store a new prospective step using the new
       * direction vector
       */
      step.propagate(function, stepLength, direction);
      observer(step.parameters.current);
    }

    // Copy the optimal parameters back into the in/out argument
    initialParameters = std::move(step.parameters.current);
    return iteration;
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
  struct CollectiveRingBuffer {
    // Ring buffered differences of parameters and gradients between updates
    Eigen::Matrix<FloatType, Eigen::Dynamic, ringBufferSize> dg;
    Eigen::Matrix<FloatType, Eigen::Dynamic, ringBufferSize> dx;
    Eigen::Matrix<FloatType, ringBufferSize, 1> dxDotDg;
    unsigned count = 0, offset = 0;

    CollectiveRingBuffer(const unsigned nParams)
      : dg(nParams, ringBufferSize),
        dx(nParams, ringBufferSize)
    {}

    void addInformation(
      const detail::EigenUpdateBuffer<VectorType>& parameterBuffer,
      const detail::EigenUpdateBuffer<VectorType>& gradientBuffer
    ) {
      if(count < ringBufferSize) {
        dg.col(count) = gradientBuffer.current - gradientBuffer.old;
        dx.col(count) = parameterBuffer.current - parameterBuffer.old;
        dxDotDg(count) = dx.col(count).dot(dg.col(count));
        ++count;
      } else {
        // Rotate columns without copying (introduce modulo offset)
        const unsigned columnOffset = (count + offset) % ringBufferSize;
        dg.col(columnOffset) = gradientBuffer.current - gradientBuffer.old;
        dx.col(columnOffset) = parameterBuffer.current - parameterBuffer.old;
        dxDotDg(columnOffset) = dx.col(columnOffset).dot(dg.col(columnOffset));
        offset = (offset + 1) % ringBufferSize;
      }
    }

    void improveDirection(Eigen::Ref<VectorType> direction) {
      VectorType alpha(count);
      auto newestToOldest = [&](unsigned i) {
        FloatType dxDotdg = dxDotDg(i);
        if(std::fabs(dxDotdg) < FloatType {1.0e-6}) {
          alpha[i] = dx.col(i).dot(direction) / FloatType {1.0e-6};
        } else {
          alpha[i] = dx.col(i).dot(direction) / dxDotdg;
        }

        direction.noalias() -= alpha[i] * dg.col(i);
      };

      // Newest to oldest
      const int newestOffset = (count + offset - 1) % ringBufferSize;
      for(int i = newestOffset; i > - 1; --i) {
        newestToOldest(i);
      }
      for(
        int i = std::min(count - 1, ringBufferSize - 1);
        i > newestOffset;
        --i
      ) {
        newestToOldest(i);
      }

      direction *= dxDotDg(newestOffset) / dg.col(newestOffset).squaredNorm();

      auto oldestToNewest = [&](unsigned i) {
        FloatType dxDotdg = dxDotDg(i);
        FloatType beta = dg.col(i).dot(direction);

        if(std::fabs(dxDotdg) < FloatType {1.0e-6}) {
          beta /= FloatType {1.0e-6};
        } else {
          beta /= dxDotDg(i);
        }

        direction.noalias() += (alpha[i] - beta) * dx.col(i);
      };

      // Oldest to newest
      const unsigned oldestOffset = (count + offset) % ringBufferSize;
      for(unsigned i = oldestOffset; i < count; ++i) {
        oldestToNewest(i);
      }
      for(unsigned i = 0; i < oldestOffset; ++i) {
        oldestToNewest(i);
      }
    }

    void updateAndGenerateNewDirection(
      Eigen::Ref<VectorType> direction,
      const StepValues& step
    ) {
      /* L-BFGS update: Use old parameters and gradients to improve steepest
       * descent direction
       */
      addInformation(step.parameters, step.gradients);
      direction = -step.gradients.current;
      improveDirection(direction);
    }
  };
};

} // namespace temple

#endif
