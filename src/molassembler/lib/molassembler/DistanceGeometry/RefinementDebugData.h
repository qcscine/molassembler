#ifndef INCLUDE_MOLASSEMBLER_DG_REFINEMENT_DEBUG_DATA_H
#define INCLUDE_MOLASSEMBLER_DG_REFINEMENT_DEBUG_DATA_H

#include "dlib/matrix.h"

#include "molassembler/DistanceGeometry/DistanceGeometry.h"

#include <list>
#include <vector>


/*! @file
 *
 * Contains the declarations of a number of debug data structures for the
 * debugging of DG refinement.
 */

namespace molassembler {

namespace DistanceGeometry {

struct RefinementStepData {
  dlib::matrix<double, 0, 1> positions;
  double distanceError;
  double chiralError;
  double fourthDimError;
  dlib::matrix<double, 0, 1> gradient;
  double proportionCorrectChiralityConstraints;
  bool compress;

  RefinementStepData(
    dlib::matrix<double, 0, 1> passPositions,
    const double passDistanceError,
    const double passChiralError,
    const double passFourthDimError,
    dlib::matrix<double, 0, 1> passGradient,
    const double passProportionCorrectChiralityConstraints,
    const bool& passCompress
  ) : positions(std::move(passPositions)),
      distanceError(passDistanceError),
      chiralError(passChiralError),
      fourthDimError(passFourthDimError),
      gradient(std::move(passGradient)),
      proportionCorrectChiralityConstraints(passProportionCorrectChiralityConstraints),
      compress(passCompress)
  {}
};

struct RefinementData {
  std::list<RefinementStepData> steps;
  std::vector<ChiralityConstraint> constraints;
  double looseningFactor;
  bool isFailure;
  std::string spatialModelGraphviz;
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
