#ifndef INCLUDE_DG_DLIB_DEBUG_ADAPTORS_H
#define INCLUDE_DG_DLIB_DEBUG_ADAPTORS_H

#include "RefinementProblem.h"
#include "DistanceGeometry/RefinementDebugData.h"

/* NOTES
 * - Identical implementations as dlibAdaptors, except that they are initialized
 *   differently with references to the refinement value functor and the list of
 *   steps that they should write debug data to
 */

namespace MoleculeManip {

namespace DistanceGeometry {

namespace dlibAdaptors {

struct debugIterationOrAllChiralitiesCorrectStrategy {
/* State */
  const unsigned maxIterations = 0;
  unsigned iterations = 0;
  // Side effects
  std::list<detail::RefinementStepData>& refinementSteps;
  // Outside state
  const errfValue<false>& valueFunctor;


/* Constructors */
  debugIterationOrAllChiralitiesCorrectStrategy(
    const unsigned& maxIter,
    std::list<detail::RefinementStepData>& refinementSteps,
    errfValue<false>& valueFunctor
  ) : maxIterations(maxIter),
      refinementSteps(refinementSteps),
      valueFunctor(valueFunctor)
  {}


/* Information */
  /*! This non-camel-case function name is necessary in order to work correctly
   * with dlib's optimizations
   */
  template<typename T>
  bool should_continue_search(
    const T& positions,
    const double& function_value __attribute__((unused)),
    const T& gradient __attribute__((unused))
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

struct debugIterationOrGradientNormStopStrategy {
/* Public access constants */
  const unsigned maxIterations;
  const double gradientNormThresholdSquared;

/* State */
  unsigned iterations = 0;
  // Side effects
  std::list<detail::RefinementStepData>& refinementSteps;
  // Outside state
  const errfValue<true>& valueFunctor;

/* Constructors */

  debugIterationOrGradientNormStopStrategy(
    const unsigned& maxIterations,
    const double& gradientNormThreshold,
    std::list<detail::RefinementStepData>& refinementSteps,
    errfValue<true>& valueFunctor
  ) : maxIterations(maxIterations),
      gradientNormThresholdSquared(gradientNormThreshold * gradientNormThreshold),
      refinementSteps(refinementSteps),
      valueFunctor(valueFunctor)
  {}

  template<typename T>
  bool should_continue_search(
    const T& positions,
    const double& function_value __attribute__((unused)),
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

} // namespace MoleculeManip

#endif
