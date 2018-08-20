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
    const dlib::matrix<double, 0, 1>& positions,
    const double& distanceError,
    const double& chiralError,
    const double& fourthDimError,
    const dlib::matrix<double, 0, 1>& gradient,
    const double& proportionCorrectChiralityConstraints,
    const bool& compress
  ) : positions(positions),
      distanceError(distanceError),
      chiralError(chiralError),
      fourthDimError(fourthDimError),
      gradient(gradient),
      proportionCorrectChiralityConstraints(proportionCorrectChiralityConstraints),
      compress(compress)
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
