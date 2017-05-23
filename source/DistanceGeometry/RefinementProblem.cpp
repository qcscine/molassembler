#include "DistanceGeometry/RefinementProblem.h"

namespace MoleculeManip {

namespace DistanceGeometry {

namespace errfDetail {
  double proportionChiralityConstraintsCorrectSign(
    const std::vector<ChiralityConstraint>& chiralityConstraints,
    const Vector& positions
  ) {
    unsigned nonZeroChiralityConstraints = 0;
    unsigned incorrectNonZeroChiralityConstraints = 0;

    for(const auto& chiralityConstraint : chiralityConstraints) {
      if(std::fabs(chiralityConstraint.lower) > 1e-4) {
        nonZeroChiralityConstraints += 1;

        const auto currentVolume = errfDetail::volume(
          positions,
          chiralityConstraint.indices
        );

        if( // can this be simplified? -> sign bit XOR?
          ( currentVolume < 0 && chiralityConstraint.lower > 0)
          || (currentVolume > 0 && chiralityConstraint.lower < 0)
        ) {
          incorrectNonZeroChiralityConstraints += 1;
        }

      }
    }

    if(nonZeroChiralityConstraints == 0) {
      return 1;
    }

    return static_cast<double>(
      nonZeroChiralityConstraints - incorrectNonZeroChiralityConstraints
    ) / nonZeroChiralityConstraints;
  }

  bool finalStructureAcceptable(
    const DistanceBoundsMatrix& bounds,
    const std::vector<ChiralityConstraint> chiralityConstraints,
    const Vector& positions
  ) {
    const double deviationThreshold = 0.5;

    // Check distance bound deviations
    for(unsigned i = 0; i < bounds.N; i++) {
      for(unsigned j = i + 1; j < bounds.N; j++) {
        double ijDistance = dlib::length(
          errfDetail::getPos(positions, j) - errfDetail::getPos(positions, i)
        );

        if(
          ijDistance - bounds.upperBound(i, j) > deviationThreshold
          || bounds.lowerBound(i, j) - ijDistance > deviationThreshold
        ) {
          return false;
        }
      }
    }

    // Check chiral bound deviations
    for(const auto& constraint : chiralityConstraints) {
      double volume = errfDetail::volume(positions, constraint.indices);

      if(
        volume - constraint.upper > deviationThreshold
        || constraint.lower - volume > deviationThreshold
      ) {
        return false;
      }
    }

    return true;
  }

} // namespace errfDetail

} // namespace DistanceGeometry

} // namespace MoleculeManip

