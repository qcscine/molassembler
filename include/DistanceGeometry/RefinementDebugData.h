#ifndef INCLUDE_DG_REFINEMENT_DEBUG_DATA_H
#define INCLUDE_DG_REFINEMENT_DEBUG_DATA_H

#include <list>
#include <vector>
#include "DistanceGeometry/DistanceGeometry.h"

#include "dlib/matrix.h"

namespace MoleculeManip {

namespace DistanceGeometry {

namespace detail {

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
};

struct DGDebugData {
  unsigned failures = 0;
  std::list<RefinementData> refinements;
};

} // namespace detail

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
