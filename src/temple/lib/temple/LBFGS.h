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
    old = std::move(current);
    current.resize(old.size());
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
    const unsigned nParams = initialParameters.size();

    detail::EigenUpdateBuffer<VectorType> parameters;
    parameters.old = std::move(initialParameters);
    detail::FloatUpdateBuffer<FloatType> values;
    detail::EigenUpdateBuffer<VectorType> gradients;
    gradients.old.resize(nParams);
    gradients.current.resize(nParams);
    CollectiveRingBuffer ringBuffer {nParams};

    /* Steepest descent step: */
    // Calculate initial function value and gradients
    function(parameters.old, values.old, gradients.old);
    observer(parameters.old);
    // Generate direction
    VectorType stepVector = FloatType {-0.1} * stepLength * gradients.old;
    // Create new parameters
    parameters.current.noalias() = parameters.old + stepLength * stepVector;
    // Re-evaluate the function, assigning current value and gradient
    function(parameters.current, values.current, gradients.current);
    observer(parameters.current);

    unsigned iteration;
    for(
      iteration = 1;
      (
        !check.checkMaxIterations(iteration)
        && !check.checkConvergence(parameters.current, values.current, gradients.current)
      );
      ++iteration
    ) {
      /* The chosen direction and step length may overshoot or undershoot,
       * and the function value might increase. We try to adapt the
       * step length in order to ensure the function value decreases and
       * step length is more or less optimal.
       */
      if(linesearch) {
        /* Armijo condition (Wolfe condition I) */
        bool armijoCondition = (
          values.current
          <= values.old + c1 * stepLength * gradients.old.dot(stepVector)
        );

        /* Curvature condition (Wolfe condition II) */
        bool curvatureCondition = (
          gradients.current.dot(stepVector)
          >= c2 * gradients.old.dot(stepVector)
        );

        /* Backtracking condition */
        bool backtrack = (values.current > values.old);

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

        /* Backtrack if we have to, otherwise accept the unoptimal, but
         * still value-improving, step
         */
        if(backtrack && (!armijoCondition || !curvatureCondition)) {
          // Retry shortened/lengthened step (no propagation needed)
          parameters.current.noalias() = parameters.old + stepLength * stepVector;
          /* Update gradients and value */
          function(parameters.current, values.current, gradients.current);
          observer(parameters.current);
          continue;
        }
      }

      /* Add a new element to the ring buffer */
      ringBuffer.addNew(parameters, gradients);

      // L-BFGS update: Use old parameters and gradients to improve direction
      stepVector = -gradients.current;
      ringBuffer.updateStepVector(stepVector);

      /* Generate new parameters and get new value and gradient */
      parameters.propagate();
      values.propagate();
      gradients.propagate();
      parameters.current.noalias() = parameters.old + stepLength * stepVector;
      function(parameters.current, values.current, gradients.current);
      observer(parameters.current);
    }

    // Copy the optimal parameters back into the in/out argument
    initialParameters = std::move(parameters.current);
    return iteration;
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

    void addNew(
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

    void updateStepVector(Eigen::Ref<VectorType> stepVector) {
      VectorType alpha(count);
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
      const unsigned oldestOffset = (count + offset) % ringBufferSize;
      for(unsigned i = oldestOffset; i < count; ++i) {
        oldestToNewest(i);
      }
      for(unsigned i = 0; i < oldestOffset; ++i) {
        oldestToNewest(i);
      }
    }
  };
};

} // namespace temple

#endif
