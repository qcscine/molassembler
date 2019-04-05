/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Adaptors for type interconversions to/from dlib
 *
 * The numerical minimization of the error function is done with the dlib
 * library. This file contains some adapters for strategy, potentials used and
 * end criteria for the minimization.
 *
 * The minimization proceeds in two stages. In the first stage, the
 * conformation is permitted to expand into a fourth spatial dimension to aid
 * in stereopermutator inversion. Afterwards, the absolute value of the fourth
 * spatial dimension is penalized and compressed out.
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_DLIB_ADAPTORS_H
#define INCLUDE_MOLASSEMBLER_DG_DLIB_ADAPTORS_H

#include "molassembler/DistanceGeometry/DlibRefinement.h"

namespace Scine {

namespace molassembler {

namespace DistanceGeometry {

namespace dlibAdaptors {

/*!
 * Dlib minimization strategy that proceeds without compression of the fourth
 * dimension until all chiral centers have inverted to the correct orientation
 */
class IterationOrAllChiralitiesCorrectStrategy {
private:
  const ErrorFunctionValue& _valueFunctor;

public:
/* State */
  std::reference_wrapper<unsigned> iterationsRef;
  const unsigned maxIterations = 0;

/* Constructors */
  // Without max iteration limit
  explicit IterationOrAllChiralitiesCorrectStrategy(
    unsigned& iterationCounter,
    const ErrorFunctionValue& valueFunctor
  ) : _valueFunctor(valueFunctor),
      iterationsRef(iterationCounter)
  {iterationCounter = 0;}

  IterationOrAllChiralitiesCorrectStrategy(
    unsigned& iterationCounter,
    const ErrorFunctionValue& valueFunctor,
    const unsigned maxIter
  ) : _valueFunctor(valueFunctor),
      iterationsRef(iterationCounter),
      maxIterations(maxIter)
  {iterationCounter = 0;}


/* Information */
  /*!
   * @note This non-camel-case function name is necessary in order to work
   *   correctly with dlib's optimizations
   */
  template<typename T>
  bool should_continue_search(
    const T& /* positions */,
    const double /* function_value */,
    const T& /* gradient */
  ) {
    iterationsRef.get() += 1;

    if(maxIterations != 0 && iterationsRef.get() > maxIterations) {
      return false;
    }

    return _valueFunctor.proportionChiralConstraintsCorrectSign < 1;
  }

};

/*!
 * Dlib minimization strategy that proceeds with compression of the fourth
 * dimension until all the overall gradient length is below a particular
 * threshold
 */
struct IterationOrGradientNormStopStrategy {
/* Public access constants */
  const unsigned maxIterations;
  const double gradientNormThresholdSquared;

/* State */
  std::reference_wrapper<unsigned> iterationsRef;

/* Constructors */

  IterationOrGradientNormStopStrategy(
    unsigned& iterationCounter,
    const unsigned passMaxIterations,
    const double gradientNormThreshold
  ) : maxIterations(passMaxIterations),
      gradientNormThresholdSquared(gradientNormThreshold * gradientNormThreshold),
      iterationsRef(iterationCounter)
  {iterationCounter = 0;}

  template<typename T>
  bool should_continue_search(
    const T& /* positions */,
    const double /* function_value */,
    const T& gradient
  ) {
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
