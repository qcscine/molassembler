#ifndef INCLUDE_MOLASSEMBLER_DG_DLIB_ADAPTORS_H
#define INCLUDE_MOLASSEMBLER_DG_DLIB_ADAPTORS_H

#include "molassembler/DistanceGeometry/RefinementProblem.h"
#include "molassembler/DistanceGeometry/RefinementDebugData.h"

/*! @file
 *
 * The numerical minimization of the error function is done with the dlib
 * library. This file contains some adapters for strategy, potentials used and
 * end criteria for the minimization.
 *
 * The minimization proceeds in two stages. In the first stage, the
 * conformation is permitted to expand into a fourth spatial dimension to aid
 * in stereocenter inversion. Afterwards, the absolute value of the fourth
 * spatial dimension is penalized and compressed out.
 */

namespace molassembler {

namespace DistanceGeometry {

namespace dlibAdaptors {

/*!
 * Dlib minimization strategy that proceeds without compression of the fourth
 * dimension until all chiral centers have inverted to the correct orientation
 */
class iterationOrAllChiralitiesCorrectStrategy {
private:
  const std::vector<ChiralityConstraint>& _constraints;

public:
/* State */
  const unsigned maxIterations = 0;
  unsigned iterations = 0;

/* Constructors */
  // Without max iteration limit
  explicit iterationOrAllChiralitiesCorrectStrategy(
    const std::vector<ChiralityConstraint>& constraints
  ) : _constraints(constraints)
  {}

  iterationOrAllChiralitiesCorrectStrategy(
    const std::vector<ChiralityConstraint>& constraints,
    const unsigned maxIter
  ) : _constraints(constraints),
      maxIterations(maxIter)
  {}


/* Information */
  /*! This non-camel-case function name is necessary in order to work correctly
   * with dlib's optimizations
   */
  template<typename T>
  bool should_continue_search(
    const T& positions,
    const double /* function_value */,
    const T& /* gradient */
  ) {
    iterations += 1;

    if(maxIterations != 0 && iterations > maxIterations) {
      return false;
    }

    if(
      errfDetail::proportionChiralityConstraintsCorrectSign(
        _constraints,
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
 * threshold
 */
struct iterationOrGradientNormStopStrategy {
/* Public access constants */
  const unsigned maxIterations;
  const double gradientNormThresholdSquared;

/* State */
  unsigned iterations = 0;

/* Constructors */

  iterationOrGradientNormStopStrategy(
    const unsigned passMaxIterations,
    const double gradientNormThreshold
  ) : maxIterations(passMaxIterations),
      gradientNormThresholdSquared(gradientNormThreshold * gradientNormThreshold)
  {}

  template<typename T>
  bool should_continue_search(
    const T& /* positions */,
    const double /* function_value */,
    const T& gradient
  ) {
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
