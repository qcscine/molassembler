// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/DistanceGeometry/RefinementProblem.h"

#include "temple/Stringify.h"

#include "molassembler/Log.h"

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

template<class Visitor>
auto visitUnfulfilledConstraints(
  const DistanceBoundsMatrix& bounds,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  const std::vector<DihedralConstraint>& dihedralConstraints,
  const Vector& positions,
  Visitor visitor
) {
  // Check distance bounds deviations
  for(unsigned i = 0; i < bounds.N(); i++) {
    for(unsigned j = i + 1; j < bounds.N(); j++) {
      double ijDistance = dlib::length(
        errfDetail::getPos(positions, j) - errfDetail::getPos(positions, i)
      );

      if(
        ijDistance - bounds.upperBound(i, j) > visitor.deviationThreshold
        || bounds.lowerBound(i, j) - ijDistance > visitor.deviationThreshold
      ) {
        visitor.distanceOverThreshold(i, j, ijDistance);
        if(visitor.earlyExit) {
          return visitor.value;
        }
      }
    }
  }

  // Check chiral bound deviations
  for(const auto& constraint : chiralityConstraints) {
    double volume = errfDetail::adjustedSignedVolume(positions, constraint.sites);

    if(
      volume - constraint.upper > visitor.deviationThreshold
      || constraint.lower - volume > visitor.deviationThreshold
    ) {
      visitor.chiralOverThreshold(constraint, volume);
      if(visitor.earlyExit) {
        return visitor.value;
      }
    }
  }

  // Check dihedral constraints
  for(const auto& constraint : dihedralConstraints) {
    double phi = errfDetail::dihedralAngle(positions, constraint.sites);

    double constraintSumHalved = (constraint.lower + constraint.upper) / 2;

    if(phi < constraintSumHalved - M_PI) {
      phi += 2 * M_PI;
    } else if(phi > constraintSumHalved + M_PI) {
      phi -= 2 * M_PI;
    }

    double term = std::fabs(phi - constraintSumHalved) - (constraint.upper - constraint.lower) / 2;

    if(term > visitor.deviationThreshold) {
      visitor.dihedralOverThreshold(constraint, term);
      if(visitor.earlyExit) {
        return visitor.value;
      }
    }
  }

  return visitor.value;
}

bool finalStructureAcceptable(
  const DistanceBoundsMatrix& bounds,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  const std::vector<DihedralConstraint>& dihedralConstraints,
  const Vector& positions
) {
  struct FinalStructureAcceptableVisitor {
    const double deviationThreshold = 0.5;
    bool earlyExit = false;
    bool value = true;

    void distanceOverThreshold(AtomIndex /* i */, AtomIndex /* j */, double /* distance */) {
      earlyExit = true;
      value = false;
    }

    void chiralOverThreshold(const ChiralityConstraint& /* chiral */, double /* volume */) {
      earlyExit = true;
      value = false;
    }

    void dihedralOverThreshold(const DihedralConstraint& /* dihedral */, double /* term */) {
      earlyExit = true;
      value = false;
    }
  };

  return visitUnfulfilledConstraints(
    bounds,
    chiralityConstraints,
    dihedralConstraints,
    positions,
    FinalStructureAcceptableVisitor {}
  );
}

void explainAcceptanceFailure(
  const DistanceBoundsMatrix& bounds,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  const std::vector<DihedralConstraint>& dihedralConstraints,
  const Vector& positions
) {
  struct AcceptanceFailureExplainer {
    const double deviationThreshold = 0.5;
    bool earlyExit = false;
    bool value = true;
    std::ostream& log = Log::log(Log::Particulars::DGStructureAcceptanceFailures);
    const DistanceBoundsMatrix& distanceBounds;

    AcceptanceFailureExplainer(const DistanceBoundsMatrix& passBounds)
      : distanceBounds(passBounds) {}

    void distanceOverThreshold(AtomIndex i, AtomIndex j, double distance) {
      log << "Distance constraint " << i << " - " << j << " : ["
        << distanceBounds.lowerBound(i, j) << ", " << distanceBounds.upperBound(i, j)
        << "] deviation over threshold, is : " << distance << "\n";
    }

    void chiralOverThreshold(const ChiralityConstraint& chiral, double volume) {
      log << "Chiral constraint " << temple::stringify(chiral.sites) << " : ["
        << chiral.lower << ", " << chiral.upper
        << "] deviation over threshold, is : " << volume << "\n";
    }

    void dihedralOverThreshold(const DihedralConstraint& dihedral, double term) {
      log << "Dihedral constraint " << temple::stringify(dihedral.sites)
        << " : [" << dihedral.lower << ", " << dihedral.upper
        << "], term is : " << term << "\n";
    }
  };

  visitUnfulfilledConstraints(
    bounds,
    chiralityConstraints,
    dihedralConstraints,
    positions,
    AcceptanceFailureExplainer {bounds}
  );
}

void explainFinalContributions(
  const DistanceBoundsMatrix& bounds,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  const std::vector<DihedralConstraint>& dihedralConstraints,
  const Vector& positions
) {
  struct FinalContributionsExplainer {
    const double deviationThreshold = 0;
    bool earlyExit = false;
    bool value = true;
    std::ostream& log = Log::log(Log::Particulars::DGFinalErrorContributions);
    const DistanceBoundsMatrix& distanceBounds;

    FinalContributionsExplainer(const DistanceBoundsMatrix& passBounds)
      : distanceBounds(passBounds) {}

    void distanceOverThreshold(AtomIndex i, AtomIndex j, double distance) {
      log << "Distance constraint " << i << " - " << j << " : ["
        << distanceBounds.lowerBound(i, j) << ", " << distanceBounds.upperBound(i, j)
        << "] is unsatisfied: " << distance << "\n";
    }

    void chiralOverThreshold(const ChiralityConstraint& chiral, double volume) {
      log << "Chiral constraint " << temple::stringify(chiral.sites) << " : ["
        << chiral.lower << ", " << chiral.upper
        << "] is unsatisfied: " << volume << "\n";
    }

    void dihedralOverThreshold(const DihedralConstraint& dihedral, double term) {
      log << "Dihedral constraint " << temple::stringify(dihedral.sites)
        << " : [" << dihedral.lower << ", " << dihedral.upper
        << "] is unsatisfied, term is: " << term << "\n";
    }
  };

  visitUnfulfilledConstraints(
    bounds,
    chiralityConstraints,
    dihedralConstraints,
    positions,
    FinalContributionsExplainer {bounds}
  );
}

} // namespace errfDetail

} // namespace DistanceGeometry

} // namespace molassembler
