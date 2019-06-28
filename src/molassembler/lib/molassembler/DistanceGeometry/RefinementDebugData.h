/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Data struct for DG refinement step data
 *
 * Contains the declarations of a number of debug data structures for the
 * debugging of DG refinement.
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_REFINEMENT_DEBUG_DATA_H
#define INCLUDE_MOLASSEMBLER_DG_REFINEMENT_DEBUG_DATA_H

#include "molassembler/DistanceGeometry/DistanceGeometry.h"

#include <Eigen/Core>

#include <list>
#include <vector>

namespace Scine {

namespace molassembler {

namespace DistanceGeometry {

struct RefinementStepData {
  Eigen::VectorXd positions;
  double distanceError;
  double chiralError;
  double dihedralError;
  double fourthDimError;
  Eigen::VectorXd gradient;
  double proportionCorrectChiralConstraints;
  bool compress;

  RefinementStepData(
    const Eigen::VectorXd& passPositions,
    const double passDistanceError,
    const double passChiralError,
    const double passDihedralError,
    const double passFourthDimError,
    const Eigen::VectorXd& passGradient,
    const double passProportionCorrectChiralConstraints,
    const bool passCompress
  ) : positions(passPositions),
      distanceError(passDistanceError),
      chiralError(passChiralError),
      dihedralError(passDihedralError),
      fourthDimError(passFourthDimError),
      gradient(passGradient),
      proportionCorrectChiralConstraints(passProportionCorrectChiralConstraints),
      compress(passCompress)
  {}
};

struct RefinementData {
  std::list<RefinementStepData> steps;
  std::vector<ChiralConstraint> constraints;
  double looseningFactor;
  bool isFailure;
  std::string spatialModelGraphviz;
};

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine

#endif
