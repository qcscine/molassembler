#ifndef INCLUDE_DG_REFINEMENT_PROBLEM_H
#define INCLUDE_DG_REFINEMENT_PROBLEM_H

#include "cppoptlib/problem.h"
#include "common_typedefs.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"

// only as long as callback exists in its current form
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

/* TODO
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
/* State */
  const std::vector<ChiralityConstraint>& _constraints;
  const DistanceBoundsMatrix& _bounds;

/* Private member functions */
  //! Make an Eigen Vector3d of an atomic index.
  const Eigen::Vector3d _getEigen(const TVector& v, const AtomIndexType& index) {
    assert(v.size() > 3 * index + 2);

    // TODO Maybe I can just get a reference to a subset instead of allocating new?
    // -> YES! member function .segment<3>(startingIndex)
    Eigen::Vector3d retv;
    retv << v(3 * index), v(3 * index + 1), v(3 * index + 2);
    return retv;
  }

  inline T _square(const T& value) { 
    return value * value;
  }

  /*!
   * Reduced since the factor 6 is not multiplied here to the result of the
   * vector dot product.
   */
  T _getTetrahedronReducedVolume(
    const TVector& v,
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k,
    const AtomIndexType& l
  ) {
    auto vecL = _getEigen(v, l);

    return (
      _getEigen(v, i) - vecL
    ).dot(
      (
        _getEigen(v, j) - vecL
      ).cross(
        _getEigen(v, k) - vecL
      )
    );
  }

  //!  First term of gradient expansion
  const Eigen::Vector3d _A(
    const TVector& v,
    const AtomIndexType& i, 
    const AtomIndexType& j
  ) {
    Eigen::Vector3d retv = _getEigen(v, j) - _getEigen(v, i);

    retv *= 4 * (
      retv.squaredNorm() / _square(
        _bounds.upperBound(i, j)
      ) - 1
    ) / _square(
      _bounds.upperBound(i, j)
    );

    return retv;
  }

  //! Second term of gradient expansion
  const Eigen::Vector3d _B(
    const TVector& v,
    const AtomIndexType& i,
    const AtomIndexType& j
  ) {
    Eigen::Vector3d retv = _getEigen(v, j) - _getEigen(v, i);

    retv *= 8 * _square(
      _bounds.lowerBound(i, j)
    ) * (
      _square(
        _bounds.lowerBound(i, j)
      ) - retv.squaredNorm()
    ) / std::pow(
      _square(
        _bounds.lowerBound(i, j)
      ) + retv.squaredNorm(),
      3
    );

    return retv;
  }

  //! Third term of gradient expansion
  T _C(
    const TVector& v,
    const AtomIndexType& i, 
    const AtomIndexType& j, 
    const AtomIndexType& k, 
    const AtomIndexType& l, 
    const double& target
  ) {
    return (-2) * (
      target - _getTetrahedronReducedVolume(v, i, j, k, l)
    );
  }

public:

  DGRefinementProblem(
    const std::vector<ChiralityConstraint>& constraints,
    const DistanceBoundsMatrix& bounds
  ) : 
    _constraints(constraints), 
    _bounds(bounds) 
  {}

  T distanceError(const TVector& v) {
    T error = 0, distance;
    const AtomIndexType N = v.size() / 3;

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = i + 1; j < N; j++) {

        distance = (
          _getEigen(v, j) - _getEigen(v, i)
        ).squaredNorm();

        // first term
        error += _square(
          distance / _square(
            _bounds.upperBound(i, j)
          ) - 1
        );

        // second term
        error += _square(
          2 * _square(
            _bounds.lowerBound(i, j)
          ) / (
            _square(
              _bounds.lowerBound(i, j)
            ) + distance
          ) - 1
        );
      }
    }

    return error;
  }

  T chiralError(const TVector& v) {
    T sum = 0;
    AtomIndexType a, b, c, d;
    T targetVal;

    for(const auto& chiralityConstraint: _constraints) {
      std::tie(a, b, c, d, targetVal) = chiralityConstraint;
      sum += _square( 
        targetVal - _getTetrahedronReducedVolume(v, a, b, c, d)
      );
    }

    return sum;
  }

  /*! 
   * Required for cppoptlib: Calculates value for specified coordinates.
   */
  T value(const TVector& v) {
    assert(v.size() % 3 == 0);

    return distanceError(v) + chiralError(v);
  }

  /*!
   * Required for cppoptlib: Calculates gradient for specified coordinates.
   */
  void gradient(const TVector& v, TVector& grad) {
    // this is where it gets VERY interesting
    const AtomIndexType N = v.size() / 3;
    Eigen::Vector3d localGradient;

    // loop temp variables
    AtomIndexType a, b, c;
    double target;
    for(AtomIndexType alpha = 0; alpha < N; alpha++) {
      localGradient.setZero();

      // A and B
      /* Skipping i == alpha is exceptionally important, for the reason that 
       * if i == alpha, although mathematically it would be assumed that
       * _A(i, i) == 0 and _B(i, i) == 0, nevertheless both the upper and 
       * lower bound for this situation are 0, leading to zeroes on 
       * denominators, leading to NaNs!
       *
       * Two versions of the loop avoidance are below, need benchmarking (TODO)
       */
      /* Version A 
      for(AtomIndexType i = 0; i < N; i++) {
        if(i == alpha) continue; 
        localGradient += _A(v, i, alpha) + _B(v, alpha, i);
      } */

      /* Version B */
      for(AtomIndexType i = 0; i < alpha; i++) {
        localGradient += _A(v, i, alpha) + _B(v, alpha, i);
      }
      for(AtomIndexType i = alpha + 1; i < N; i++) {
        localGradient += _A(v, i, alpha) + _B(v, alpha, i);
      }

      /* C 
       * (chirality constraints)
       * TODO below: why and how?
       */
      for(const auto& chiralityConstraint: _constraints) {
        if(std::get<0>(chiralityConstraint) == alpha) {
          std::tie(std::ignore, a, b, c, target) = chiralityConstraint;
        } else if(std::get<1>(chiralityConstraint) == alpha) {
          std::tie(a, std::ignore, b, c, target) = chiralityConstraint;
        } else if(std::get<2>(chiralityConstraint) == alpha) {
          std::tie(a, b, std::ignore, c, target) = chiralityConstraint;
        } else if(std::get<3>(chiralityConstraint) == alpha) {
          std::tie(a, b, c, std::ignore, target) = chiralityConstraint;
        }

        localGradient += _C(v, alpha, a, b, c, target) * (
          _getEigen(v, a) - _getEigen(v, c)
        ).cross(
          _getEigen(v, b) - _getEigen(v, c)
        );
      }

      // Assign to gradient
      grad(3 * alpha) = localGradient(0);
      grad(3 * alpha + 1) = localGradient(1);
      grad(3 * alpha + 2) = localGradient(2);
    }
  }

  bool callback(
    const cppoptlib::Criteria<T>& state,
    const TVector& x
  ) {
    std::cout << "(" << std::setw(2) << state.iterations << ")"
              << " ||dx|| = " << std::fixed << std::setw(8) << std::setprecision(4) << state.gradNorm
              << " ||x|| = "  << std::setw(6) << x.norm()
              << " f(x) = "   << std::setw(8) << value(x)
              << " x = [" << std::setprecision(8) << x.transpose() << "]" << std::endl;
    return true;
  }

};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
