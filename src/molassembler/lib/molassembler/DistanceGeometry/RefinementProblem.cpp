#include "DistanceGeometry/RefinementProblem.h"
#include "Log.h"
#include "temple/Stringify.h"

namespace molassembler {

namespace DistanceGeometry {

namespace errfDetail {
  double proportionChiralityConstraintsCorrectSign(
    const std::vector<ChiralityConstraint>& chiralityConstraints,
    const Vector& positions
  ) {
    unsigned nonZeroChiralityConstraints = 0;
    unsigned incorrectNonZeroChiralityConstraints = 0;

    for(const auto& chiralityConstraint : chiralityConstraints) {
      /* Make sure that chirality constraints meant to cause coplanarity of
       * four indices aren't counted here - Their volume bounds are
       * usually so tight that they might be slightly outside in any given
       * refinement step. If they are outside their target values by a large
       * amount, that is caught by finalStructureAcceptable anyway.
       */
      if(std::fabs(chiralityConstraint.lower + chiralityConstraint.upper) > 1e-4) {
        nonZeroChiralityConstraints += 1;

        const auto currentVolume = errfDetail::adjustedSignedVolume(
          positions,
          chiralityConstraint.sites
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
    for(unsigned i = 0; i < bounds.N(); i++) {
      for(unsigned j = i + 1; j < bounds.N(); j++) {
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
      double volume = errfDetail::adjustedSignedVolume(positions, constraint.sites);

      if(
        volume - constraint.upper > deviationThreshold
        || constraint.lower - volume > deviationThreshold
      ) {
        return false;
      }
    }

    return true;
  }

  void explainAcceptanceFailure(
    const DistanceBoundsMatrix& bounds,
    const std::vector<ChiralityConstraint> chiralityConstraints,
    const Vector& positions
  ) {
    auto& log = Log::log(Log::Particulars::DGStructureAcceptanceFailures);
    const double deviationThreshold = 0.5;

    // Check distance bound deviations
    for(unsigned i = 0; i < bounds.N(); i++) {
      for(unsigned j = i + 1; j < bounds.N(); j++) {
        double ijDistance = dlib::length(
          errfDetail::getPos(positions, j) - errfDetail::getPos(positions, i)
        );

        if(
          ijDistance - bounds.upperBound(i, j) > deviationThreshold
          || bounds.lowerBound(i, j) - ijDistance > deviationThreshold
        ) {
          log << "Distance constraint " << i << " - " << j << " : ["
            << bounds.lowerBound(i, j) << ", " << bounds.upperBound(i, j)
            << "] deviation over threshold, is : " << ijDistance << "\n";
        }
      }
    }

    // Check chiral bound deviations
    for(const auto& constraint : chiralityConstraints) {
      double volume = errfDetail::adjustedSignedVolume(positions, constraint.sites);

      if(
        volume - constraint.upper > deviationThreshold
        || constraint.lower - volume > deviationThreshold
      ) {
        log << "Chiral constraint " << temple::stringify(constraint.sites) << " : ["
          << constraint.lower << ", " << constraint.upper
          << "] deviation over threshold, is : " << volume << "\n";
      } else if(
        volume - constraint.upper > 0
        || constraint.lower - volume > 0
      ) {
        log << "Unsatisfied chiral constraint (within threshold for acceptance) "
          << temple::stringify(constraint.sites) << ": ["
          << constraint.lower << ", " << constraint.upper
          << "], is : " << volume << "\n";
      }
    }
  }

} // namespace errfDetail

} // namespace DistanceGeometry

} // namespace molassembler
