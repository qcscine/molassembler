#ifndef INCLUDE_MOLASSEMBLER_DG_DLIB_DEBUG_ADAPTORS_H
#define INCLUDE_MOLASSEMBLER_DG_DLIB_DEBUG_ADAPTORS_H

#include "molassembler/DistanceGeometry/RefinementDebugData.h"
#include "molassembler/DistanceGeometry/RefinementProblem.h"

/*! @file
 *
 * @brief dlib DG refinement minimization strategies for debugging
 *
 * The numerical minimization of the error function is done with the dlib
 * library. This file contains some adapters for strategy, potentials used and
 * end criteria for the minimization.
 *
 * The minimization proceeds in two stages. In the first stage, the
 * conformation is permitted to expand into a fourth spatial dimension to aid
 * in stereocenter inversion. Afterwards, the absolute value of the fourth
 * spatial dimension is penalized and compressed out.
 *
 * This set of minimization strategies is intended for debugging. They create,
 * as a side effect, lits of refinement step data.
 */

/* NOTES
 * - Identical implementations as dlibAdaptors, except that they are initialized
 *   differently with references to the refinement value functor and the list of
 *   steps that they should write debug data to
 */

namespace molassembler {

namespace DistanceGeometry {

namespace dlibAdaptors {

/*!
 * Dlib minimization strategy that proceeds without compression of the fourth
 * dimension until all chiral centers have inverted to the correct orientation.
 *
 * Creates a list of refinement step data via side-effect.
 */
struct debugIterationOrAllChiralitiesCorrectStrategy {
/* State */
  const unsigned maxIterations = 0;
  unsigned iterations = 0;
  // Side effects
  std::list<RefinementStepData>& refinementSteps;
  // Outside state
  const errfValue<false>& valueFunctor;


/* Constructors */
  debugIterationOrAllChiralitiesCorrectStrategy(
    const unsigned maxIter,
    std::list<RefinementStepData>& passRefinementSteps,
    errfValue<false>& passValueFunctor
  ) : maxIterations(maxIter),
      refinementSteps(passRefinementSteps),
      valueFunctor(passValueFunctor)
  {}


/* Information */
  /*! This non-camel-case function name is necessary in order to work correctly
   * with dlib's optimizations
   */
  template<typename T>
  bool should_continue_search(
    const T& positions,
    const double /* function_value */,
    const T& gradient
  ) {
    refinementSteps.emplace_back(
      positions,
      valueFunctor.distanceError(positions),
      valueFunctor.chiralError(positions),
      valueFunctor.extraDimensionError(positions),
      gradient,
      errfDetail::proportionChiralityConstraintsCorrectSign(
        valueFunctor.chiralityConstraints,
        positions
      ),
      false
    );

    iterations += 1;

    if(maxIterations != 0 && iterations > maxIterations) {
      return false;
    }

    if(
      errfDetail::proportionChiralityConstraintsCorrectSign(
        valueFunctor.chiralityConstraints,
        positions
      ) >= 1 // just in case >
    ) {
      return false;
    }

    return true;
  }
};

/*!
 * Dlib minimization strategy that proceeds with compression of the fourth
 * dimension until all the overall gradient length is below a particular
 * threshold.
 *
 * Creates a list of refinement step data as side effect.
 */
struct debugIterationOrGradientNormStopStrategy {
/* Public access constants */
  const unsigned maxIterations;
  const double gradientNormThresholdSquared;

/* State */
  unsigned iterations = 0;
  // Side effects
  std::list<RefinementStepData>& refinementSteps;
  // Outside state
  const errfValue<true>& valueFunctor;

/* Constructors */

  debugIterationOrGradientNormStopStrategy(
    const unsigned passMaxIterations,
    const double gradientNormThreshold,
    std::list<RefinementStepData>& passRefinementSteps,
    const errfValue<true>& passValueFunctor
  ) : maxIterations(passMaxIterations),
      gradientNormThresholdSquared(gradientNormThreshold * gradientNormThreshold),
      refinementSteps(passRefinementSteps),
      valueFunctor(passValueFunctor)
  {}

  template<typename T>
  bool should_continue_search(
    const T& positions,
    const double /* function_value */,
    const T& gradient
  ) {
    refinementSteps.emplace_back(
      positions,
      valueFunctor.distanceError(positions),
      valueFunctor.chiralError(positions),
      valueFunctor.extraDimensionError(positions),
      gradient,
      errfDetail::proportionChiralityConstraintsCorrectSign(
        valueFunctor.chiralityConstraints,
        positions
      ),
      true
    );

    iterations += 1;

    if(maxIterations != 0 && iterations > maxIterations) {
      return false;
    }

    /* Optimization:
     * Instead of comparing
     *   dlib::length(gradient) < gradientNormThreshold
     * = std::sqrt(dlib::length_squared(gradient)) < gradientNormThreshold
     *
     * we square both sides, and since both quantities are > 0,
     *  dlib::length_squared(gradient) < gradientNormThreshold²
     *
     * This saves a std::sqrt every iteration
     */
    if(dlib::length_squared(gradient) < gradientNormThresholdSquared) {
      return false;
    }

    return true;
  }
};

} // namespace dlibAdaptors

} // namespace DistanceGeometry

} // namespace molassembler

#endif
