#ifndef INCLUDE_DG_REFINEMENT_PROBLEM_H
#define INCLUDE_DG_REFINEMENT_PROBLEM_H

#include "cppoptlib/problem.h"
#include "common_typedefs.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "Log.h"

// only as long as callback exists in its current form
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

/* TODO
 * - Correctness tests
 * - Optimize 
 *   - remove all the bound squaring be pre-squaring the bounds
 *   - ... ?
 */

namespace MoleculeManip {

namespace DistanceGeometry {

/*! Although it's templated, this is only instantiable as a double! T must match
 * the type of the target value in a ChiralityConstraint -> See the std::tie 
 * call when chirality constraints are considered
 */
template<typename T>
class DGRefinementProblem : public cppoptlib::Problem<T> {
public:
/* Typedefs */
  using typename cppoptlib::Problem<T>::TVector;
  using typename cppoptlib::Problem<T>::THessian;

private:
/* Private member functions */
  //! Make an Eigen Vector4d of an atomic index.
  template<int vectorSize = 4>
  auto _getEigen(const TVector& v, const AtomIndexType& index) const {
    assert(v.size() > 4 * index + 3);

    /* Return a fixed-size const reference to a part of the vector
     *
     * Interesting tidbit to the syntax below:
     * If you write: return v.segment<3>(3 * index);
     *             member function -^^^- integer 3
     *                               |
     *                          operator <
     * 
     * so to disambiguate, you write template before the name of the member
     * function.
     */
    return v.template segment<vectorSize>(4 * index);
  }

  inline T _square(const T& value) const { 
    return value * value;
  }

  /*!
   * Reduced since the factor 6 is not multiplied here to the result of the
   * vector dot product.
   */
  T _getTetrahedronReducedVolume(
    const TVector& v,
    const std::array<AtomIndexType, 4>& indices
  ) const {
    auto vecL = _getEigen<3>(v, indices[3]);

    return (
      _getEigen<3>(v, indices[0]) - vecL
    ).dot(
      (
        _getEigen<3>(v, indices[1]) - vecL
      ).cross(
        _getEigen<3>(v, indices[2]) - vecL
      )
    );
  }

  //!  First term of gradient expansion
  Eigen::Vector4d _A(
    const TVector& v,
    const AtomIndexType& i, 
    const AtomIndexType& j
  ) const {
    Eigen::Vector4d retv = _getEigen(v, j) - _getEigen(v, i);

    retv *= 4 * (
      retv.squaredNorm() / _square(
        bounds.upperBound(i, j)
      ) - 1
    ) / _square(
      bounds.upperBound(i, j)
    );

    return retv;
  }

  //! Second term of gradient expansion
  Eigen::Vector4d _B(
    const TVector& v,
    const AtomIndexType& i,
    const AtomIndexType& j
  ) const {
    Eigen::Vector4d retv = _getEigen(v, j) - _getEigen(v, i);

    retv *= 8 * _square(
      bounds.lowerBound(i, j)
    ) * (
      _square(
        bounds.lowerBound(i, j)
      ) - retv.squaredNorm()
    ) / std::pow(
      _square(
        bounds.lowerBound(i, j)
      ) + retv.squaredNorm(),
      3
    );

    return retv;
  }

  //! Third term of gradient expansion
  T _C(
    const TVector& v,
    const std::array<AtomIndexType, 4>& indices,
    const double& target
  ) const {
    return (-2) * (
      target - _getTetrahedronReducedVolume(v, indices)
    );
  }

public:
/* State */
  const std::vector<ChiralityConstraint>& constraints;
  const DistanceBoundsMatrix& bounds;
  mutable bool compress = false;

/* Constructors */
  DGRefinementProblem(
    const std::vector<ChiralityConstraint>& constraints,
    const DistanceBoundsMatrix& bounds
  ) : 
    constraints(constraints), 
    bounds(bounds) 
  {}

/* Information */
  T distanceError(const TVector& v) const {
    T error = 0, distance;
    const AtomIndexType N = v.size() / 4;

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = i + 1; j < N; j++) {

        distance = (
          _getEigen(v, j) - _getEigen(v, i)
        ).squaredNorm();

        // first term
        error += _square(
          distance / _square(
            bounds.upperBound(i, j)
          ) - 1
        );

        // second term
        error += _square(
          2 * _square(
            bounds.lowerBound(i, j)
          ) / (
            _square(
              bounds.lowerBound(i, j)
            ) + distance
          ) - 1
        );
      }
    }

    return error;
  }

  void invertY(TVector& v) const {
    const unsigned N = v.size() / 4;

    for(unsigned i = 0; i < N; i++) {
      v[4 * i + 1] *= -1;
    }
  }

  double proportionCorrectChiralityConstraints(const TVector& v) const {
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

  T chiralError(const TVector& v) const {
    T sum = 0;

    for(const auto& chiralityConstraint: constraints) {
      sum += _square( 
        chiralityConstraint.target 
        - _getTetrahedronReducedVolume(v, chiralityConstraint.indices)
      );
    }

    return sum;
  }

  T extraDimensionError(const TVector& v) const {
    if(!compress) return 0;

    T sum = 0;

    const AtomIndexType N = v.size() / 4;
    for(unsigned i = 0; i < N; i++) {
      sum += _square(
        v[4 * i + 3]
      );
    }

    return sum;
  }

  /*! 
   * Calculates value for specified coordinates.
   */
  T value(const TVector& v) override {
    assert(v.size() % 4 == 0);

    return distanceError(v) + chiralError(v) + extraDimensionError(v);
  }

  /*!
   * Required for cppoptlib: Calculates gradient for specified coordinates.
   */
  void gradient(const TVector& v, TVector& grad) const {
    /* At the beginning of every gradient calculation, we check if all chiral
     * centers are correct, so we can switch on 4th-dimension compression
     */
    if(!compress && proportionCorrectChiralityConstraints(v) == 1.0) {
      compress = true;
    }

    const AtomIndexType N = v.size() / 4;
    Eigen::Vector4d localGradient;

    for(AtomIndexType alpha = 0; alpha < N; alpha++) {
      localGradient.setZero();

      // A and B
      /* Skipping i == alpha is exceptionally important, for the reason that 
       * if i == alpha, although mathematically it would be assumed that
       * _A(i, i) == 0 and _B(i, i) == 0, nevertheless both the upper and 
       * lower bound for this situation are 0, leading to zeroes on 
       * denominators, leading to NaNs!
       */
      for(AtomIndexType i = 0; i < alpha; i++) {
        localGradient += _A(v, i, alpha) + _B(v, alpha, i);
      }
      for(AtomIndexType i = alpha + 1; i < N; i++) {
        localGradient += _A(v, i, alpha) + _B(v, alpha, i);
      }

      /* C 
       * (chirality constraints)
       */
      for(const auto& constraint: constraints) {
        bool fallthrough = false;
        double target = constraint.target;
        std::array<AtomIndexType, 4> indices = constraint.indices;
        // TODO refactor this, there has to be a good way to write this
        // if(constraint.indices[0] == alpha) -> do nothing!
        if(constraint.indices[1] == alpha) {
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
          // I'm pretty confident here now that this is correct
/* #ifndef NDEBUG
          auto temp = _C(v, indices, target) * (
            _getEigen<3>(v, indices[1]) - _getEigen<3>(v, indices[3])
          ).cross(
            _getEigen<3>(v, indices[2]) - _getEigen<3>(v, indices[3])
          );

          localGradient.template segment<3>(0) += temp;

          if(temp.norm() > 10) {
            auto& streamRef = Log::log(Log::Particulars::DGRefinementChiralityNumericalDebugInfo);
            streamRef << "Unusually large chirality gradient contribution on alpha = "
              << alpha << ", constraint = vec{";

            for(unsigned i = 0; i < 4; i++) {
              streamRef << constraint.indices[i];
              if(i != 3) streamRef << ", ";
            }

            streamRef << "}, reordered indices = vec{";

            for(unsigned i = 0; i < 4; i++) {
              streamRef << indices[i];
              if(i != 3) streamRef << ", ";
            }

            streamRef << "}: _C = " << _C(v, indices, target) 
              << ", target = " << target 
              << ", current volume = " << _getTetrahedronReducedVolume(v, indices) 
              << std::endl;
          }
#else
#endif
*/
          localGradient.template segment<3>(0) += _C(v, indices, target) * (
            _getEigen<3>(v, indices[1]) - _getEigen<3>(v, indices[3])
          ).cross(
            _getEigen<3>(v, indices[2]) - _getEigen<3>(v, indices[3])
          );
        }
      }

      // Gradient due to fourth coordinate if compress is on
      if(compress) {
        localGradient[3] += 2 * v[4 * alpha + 3];
      }

      // Assign to gradient
      grad.template segment<4>(4 * alpha) = localGradient;
    }
  }

  bool callback(
    const cppoptlib::Criteria<T>& state,
    const TVector& x
  ) const {
#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::DGRefinementProgress) == 1) {
      TVector gradientVector(x.size());
      gradient(x, gradientVector);

      // CSV format
      Log::log(Log::Particulars::DGRefinementProgress) 
        << x.size()/4 << ","
        << state.iterations << "," 
        << std::fixed << std::setprecision(4) << state.gradNorm << "," 
        << distanceError(x) << ","
        << chiralError(x) << ","
        << extraDimensionError(x) << ","
        << std::setprecision(4) << x.transpose() << ","
        << std::setprecision(4) << gradientVector.transpose() 
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
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
