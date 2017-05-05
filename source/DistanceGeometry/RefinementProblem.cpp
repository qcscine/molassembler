#include "DistanceGeometry/RefinementProblem.h"

#include "Log.h"

// only as long as callback exists in its current form
#include <iostream>
#include <iomanip>

namespace MoleculeManip {

namespace DistanceGeometry {

RefinementProblem::RefinementProblem(
  const std::vector<ChiralityConstraint>& constraints,
  const DistanceBoundsMatrix& bounds
) : 
  constraints(constraints), 
  squaredBounds(bounds.makeSquaredBoundsMatrix()) 
{}

RefinementProblem::RefinementProblem(
  const std::vector<ChiralityConstraint>& constraints,
  const DistanceBoundsMatrix& bounds,
  const CallbackFunction& callback
) : 
  constraints(constraints), 
  squaredBounds(bounds.makeSquaredBoundsMatrix()),
  callbackClosure(callback)
{}

bool RefinementProblem::callback(
  const cppoptlib::Criteria<double>& state,
  const TVector& x
) const {
  // Call the closure if present
  if(callbackClosure) callbackClosure(state, x, *this);

#ifndef NDEBUG
  if(Log::particulars.count(Log::Particulars::DGRefinementProgress) == 1) {

    // CSV format
    Log::log(Log::Particulars::DGRefinementProgress) 
      << x.size()/4 << ","
      << state.iterations << "," 
      << std::fixed << std::setprecision(4) << state.gradNorm << "," 
      << distanceError(x) << ","
      << chiralError(x) << ","
      << extraDimensionError(x) << ","
      << std::setprecision(4) << x.transpose() 
      << std::endl; 

    // human-readable format
    /* std::cout << "(" << std::setw(2) << state.iterations << ")"
        << " ||dx|| = " << std::fixed << std::setw(8) 
        << std::setprecision(4) << state.gradNorm
        << " ||x|| = "  << std::setw(6) << x.norm()
        << " f(x) = "   << std::setw(8) << value(x)
        << " x = [" << std::setprecision(8) << x.transpose() << "]" 
        << std::endl; */
  }
#endif

  return true;
}

double RefinementProblem::chiralError(const TVector& v) const {
  double sum = 0;

  for(const auto& chiralityConstraint: constraints) {
    sum += _square( 
      chiralityConstraint.target 
      - _getTetrahedronReducedVolume(v, chiralityConstraint.indices)
    );
  }

  return sum;
}

double RefinementProblem::extraDimensionError(const TVector& v) const {
  if(!compress) return 0;

  double sum = 0;

  const AtomIndexType N = v.size() / 4;
  for(unsigned i = 0; i < N; i++) {
    sum += _square(
      v[4 * i + 3]
    );
  }

  return sum;
}


double RefinementProblem::distanceError(const TVector& v) const {
  double error = 0, distance;
  const AtomIndexType N = v.size() / 4;

  for(unsigned i = 0; i < N; i++) {
    for(unsigned j = i + 1; j < N; j++) {

      distance = (
        _getPos(v, j) - _getPos(v, i)
      ).squaredNorm();

      // first term
      error += _square(
        distance / squaredBounds.upperBound(i, j)
        - 1
      );

      // second term
      error += _square(
        2 * squaredBounds.lowerBound(i, j) / (
          squaredBounds.lowerBound(i, j) + distance
        ) 
        - 1
      );
    }
  }

  return error;
}

void RefinementProblem::alternateGradient(const TVector& v, TVector& gradient) const {
  gradient.setZero();
  if(!compress && proportionCorrectChiralityConstraints(v) == 1.0) {
    compress = true;
  }

  gradientA(v, gradient);
  gradientB(v, gradient);
  gradientC(v, gradient);
  gradientD(v, gradient);
}

/*!
 * Required for cppoptlib: Calculates gradient for specified coordinates.
 */
void RefinementProblem::gradient(const TVector& v, TVector& grad) const {
  /* At the beginning of every gradient calculation, we check if all chiral
   * centers are correct, so we can switch on 4th-dimension compression
   */
  if(!compress && proportionCorrectChiralityConstraints(v) == 1.0) {
    compress = true;
  }

  grad.setZero();

  const AtomIndexType N = v.size() / 4;

  for(AtomIndexType alpha = 0; alpha < N; alpha++) {
    // A and B
    /* Skipping i == alpha is exceptionally important, for the reason that if
     * i == alpha, although mathematically it would be assumed that
     * gradientTermA(i, i) == 0 and gradientTermB(i, i) == 0, nevertheless
     * both the upper and lower bound for this situation are 0, leading to
     * zeroes on denominators, leading to NaNs!
     */
    for(AtomIndexType i = 0; i < alpha; i++) {
      grad.template segment<4>(4 * alpha) += gradientTermA(v, i, alpha) + gradientTermB(v, alpha, i);
    }
    for(AtomIndexType i = alpha + 1; i < N; i++) {
      grad.template segment<4>(4 * alpha) += gradientTermA(v, i, alpha) + gradientTermB(v, alpha, i);
    }

    /* C 
     * (chirality constraints)
     */
    for(const auto& constraint: constraints) {
      bool fallthrough = false; 
      double target = constraint.target;
      std::array<AtomIndexType, 4> indices;
      // TODO refactor this, there has to be a good way to write this
      // -> find alpha in constraint. if not found, fallthrough. if found,
      //    copy_if to new array. if position at which found an uneven index
      //    (mod 2), target *= -1;
      //
      if(constraint.indices[0] == alpha) {
        indices = constraint.indices;
      } else if(constraint.indices[1] == alpha) {
        indices = {
          alpha,
          constraint.indices[0],
          constraint.indices[2],
          constraint.indices[3]
        };
        
        // uneven permutation
        target *= -1;
      } else if(constraint.indices[2] == alpha) {
        indices = {
          alpha,
          constraint.indices[0],
          constraint.indices[1],
          constraint.indices[3]
        };
      } else if(constraint.indices[3] == alpha) {
        indices = {
          alpha,
          constraint.indices[0],
          constraint.indices[1],
          constraint.indices[2]
        };

        // uneven permutation
        target *= -1;
      } else {
        fallthrough = true;
      }

      if(!fallthrough) {
        grad.template segment<3>(4 * alpha) += gradientTermC(v, indices, target) * (
          _getPos<3>(v, indices[1]) - _getPos<3>(v, indices[3])
        ).cross(
          _getPos<3>(v, indices[2]) - _getPos<3>(v, indices[3])
        );
      }
    }

    // Gradient due to fourth coordinate if compress is on
    if(compress) {
      grad[4 * alpha + 3] += 2 * v[4 * alpha + 3];
    }
  }
}

void RefinementProblem::gradientA(const TVector& v, TVector& gradient) const {
  const unsigned N = v.size() / 4;

  for(unsigned alpha = 0; alpha < N; alpha++) {
    // Once up to but not including alpha
    for(unsigned i = 0; i < alpha; i++) {
      const Eigen::Vector4d diff = _getPos(v, alpha) - _getPos(v, i);

      gradient.template segment<4>(4 * alpha) += diff * 4 * (
        diff.squaredNorm() / squaredBounds.upperBound(i, alpha)
        - 1
      ) / squaredBounds.upperBound(i, alpha);
    }

    // And from alpha + 1 the remaining ones
    for(unsigned i = alpha + 1; i < N; i++) {
      const Eigen::Vector4d diff = _getPos(v, alpha) - _getPos(v, i);

      gradient.template segment<4>(4 * alpha) += diff * 4 * (
        diff.squaredNorm() / squaredBounds.upperBound(i, alpha)
        - 1
      ) / squaredBounds.upperBound(i, alpha);
    }
  }
}

void RefinementProblem::gradientB(const TVector& v, TVector& gradient) const {
  const unsigned N = v.size() / 4;

  for(unsigned alpha = 0; alpha < N; alpha++) {
    // Once up to but not including alpha
    for(unsigned i = 0; i < alpha; i++) {
      const Eigen::Vector4d diff = _getPos(v, i) - _getPos(v, alpha);

      gradient.template segment<4>(4 * alpha) += (
        diff * 8 * squaredBounds.lowerBound(i, alpha) * ( 
          squaredBounds.lowerBound(i, alpha)
          - diff.squaredNorm()
        ) / std::pow(
          squaredBounds.lowerBound(i, alpha) + diff.squaredNorm(),
          3
        )
      );
    }

    // And from alpha + 1 the remaining ones
    for(unsigned i = alpha + 1; i < N; i++) {
      const Eigen::Vector4d diff = _getPos(v, i) - _getPos(v, alpha);

      gradient.template segment<4>(4 * alpha) += (
        diff * 8 * squaredBounds.lowerBound(i, alpha) * ( 
          squaredBounds.lowerBound(i, alpha)
          - diff.squaredNorm()
        ) / std::pow(
          squaredBounds.lowerBound(i, alpha) + diff.squaredNorm(),
          3
        )
      );
    }
  }
}

void RefinementProblem::gradientC(const TVector& v, TVector& gradient) const {
  if(!constraints.empty()) {
    const unsigned N = v.size() / 4;

    for(unsigned alpha = 0; alpha < N; alpha++) {
      for(const auto& constraint: constraints) {
        auto findIter = std::find(
          constraint.indices.begin(),
          constraint.indices.end(),
          alpha
        );

        if(findIter != constraint.indices.end()) {
          double target = constraint.target;
          if((findIter - constraint.indices.begin()) % 2 == 1) {
            target *= -1;
          }

          std::array<AtomIndexType, 4> indices = {alpha};
          std::copy_if(
            constraint.indices.begin(),
            constraint.indices.end(),
            indices.begin() + 1,
            [&](const AtomIndexType& value) -> bool {
              return value != alpha;
            }
          );

          gradient.template segment<3>(4 * alpha) -= 2 * (
            target - _getTetrahedronReducedVolume(v, indices)
          ) * (
            _getPos<3>(v, indices[1]) - _getPos<3>(v, indices[3])
          ).cross(
            _getPos<3>(v, indices[2]) - _getPos<3>(v, indices[3])
          );
        }
      }
    }
  }
}

void RefinementProblem::gradientD(const TVector& v, TVector& gradient) const {
  if(compress) {
    const unsigned N = v.size() / 4;

    for(unsigned alpha = 0; alpha < N; alpha++) {
      gradient[4 * alpha + 3] += 2 * v[4 * alpha + 3];
    }
  }
}

void RefinementProblem::invertY(TVector& v) const {
  const unsigned N = v.size() / 4;

  for(unsigned i = 0; i < N; i++) {
    v[4 * i + 1] *= -1;
  }
}

double RefinementProblem::proportionCorrectChiralityConstraints(const TVector& v) const {
  unsigned nonZeroChiralityConstraints = 0;
  unsigned incorrectNonZeroChiralityConstraints = 0;

  for(const auto& chiralityConstraint : constraints) {
    if(std::fabs(chiralityConstraint.target) > 1e-4) {
      nonZeroChiralityConstraints += 1;

      auto eval = _getTetrahedronReducedVolume(
        v,
        chiralityConstraint.indices
      );

      if( // can this be simplified? -> sign bit XOR?
        ( eval < 0 && chiralityConstraint.target > 0)
        || (eval > 0 && chiralityConstraint.target < 0)
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

/*! 
 * Calculates value for specified coordinates.
 */
double RefinementProblem::value(const TVector& v) {
  assert(v.size() % 4 == 0);

  return distanceError(v) + chiralError(v) + extraDimensionError(v);
}

} // namespace DistanceGeometry

} // namespace MoleculeManip
