/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief dlib DG refinement minimization strategies for debugging
 *
 * The numerical minimization of the error function is done with the dlib
 * library. This file contains some adapters for strategy, potentials used and
 * end criteria for the minimization.
 *
 * The minimization proceeds in two stages. In the first stage, the
 * conformation is permitted to expand into a fourth spatial dimension to aid
 * in stereopermutator inversion. Afterwards, the absolute value of the fourth
 * spatial dimension is penalized and compressed out.
 *
 * This set of minimization strategies is intended for debugging. They create,
 * as a side effect, lits of refinement step data.
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_DLIB_DEBUG_ADAPTORS_H
#define INCLUDE_MOLASSEMBLER_DG_DLIB_DEBUG_ADAPTORS_H

#include "molassembler/DistanceGeometry/RefinementDebugData.h"
#include "molassembler/DistanceGeometry/DlibRefinement.h"

/* NOTES
 * - Identical implementations as dlibAdaptors, except that they are initialized
 *   differently with references to the refinement value functor and the list of
 *   steps that they should write debug data to
 */

namespace Scine {

namespace molassembler {

namespace DistanceGeometry {

namespace dlibAdaptors {

/*!
 * Dlib minimization strategy that proceeds without compression of the fourth
 * dimension until all chiral centers have inverted to the correct orientation.
 *
 * Creates a list of refinement step data via side-effect.
 */
struct DebugIterationOrAllChiralitiesCorrectStrategy {
/* State */
  const unsigned maxIterations = 0;
  std::reference_wrapper<unsigned> iterationsRef;
  // Side effects
  std::list<RefinementStepData>& refinementSteps;
  // Outside state
  const ErrorFunctionValue& valueFunctor;


/* Constructors */
  DebugIterationOrAllChiralitiesCorrectStrategy(
    unsigned& iterationCounter,
    const unsigned maxIter,
    std::list<RefinementStepData>& passRefinementSteps,
    const ErrorFunctionValue& passValueFunctor
  ) : maxIterations(maxIter),
      iterationsRef(iterationCounter),
      refinementSteps(passRefinementSteps),
      valueFunctor(passValueFunctor)
  {iterationCounter = 0;}


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
      valueFunctor.dihedralError(positions),
      valueFunctor.extraDimensionError(positions),
      gradient,
      valueFunctor.proportionChiralConstraintsCorrectSign,
      valueFunctor.compressFourthDimension
    );

    iterationsRef.get() += 1;

    if(maxIterations != 0 && iterationsRef.get() > maxIterations) {
      return false;
    }

    return valueFunctor.proportionChiralConstraintsCorrectSign < 1;
  }
};

/*!
 * Dlib minimization strategy that proceeds with compression of the fourth
 * dimension until all the overall gradient length is below a particular
 * threshold.
 *
 * Creates a list of refinement step data as side effect.
 */
struct DebugIterationOrGradientNormStopStrategy {
/* Public access constants */
  const unsigned maxIterations;
  std::reference_wrapper<unsigned> iterationsRef;
  const double gradientNormThresholdSquared;
  // Side effects
  std::list<RefinementStepData>& refinementSteps;
  // Outside state
  const ErrorFunctionValue& valueFunctor;

/* Constructors */

  DebugIterationOrGradientNormStopStrategy(
    unsigned& iterationCounter,
    const unsigned passMaxIterations,
    const double gradientNormThreshold,
    std::list<RefinementStepData>& passRefinementSteps,
    const ErrorFunctionValue& passValueFunctor
  ) : maxIterations(passMaxIterations),
      iterationsRef(iterationCounter),
      gradientNormThresholdSquared(gradientNormThreshold * gradientNormThreshold),
      refinementSteps(passRefinementSteps),
      valueFunctor(passValueFunctor)
  {iterationCounter = 0;}

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
      valueFunctor.dihedralError(positions),
      valueFunctor.extraDimensionError(positions),
      gradient,
      valueFunctor.proportionChiralConstraintsCorrectSign,
      valueFunctor.compressFourthDimension
    );

    iterationsRef.get() += 1;

    if(maxIterations != 0 && iterationsRef.get() > maxIterations) {
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

} // namespace Scine

#endif
