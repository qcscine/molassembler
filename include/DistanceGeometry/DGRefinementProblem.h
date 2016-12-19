#ifndef INCLUDE_DG_REFINEMENT_PROBLEM_H
#define INCLUDE_DG_REFINEMENT_PROBLEM_H

#include "cppoptlib/problem.h"
#include "common_typedefs.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"

// only as long as callback exists in its current form
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

namespace MoleculeManip {

namespace DistanceGeometry {

template<typename T>
class DGRefinementProblem : public cppoptlib::Problem<T> {
public:
  using typename cppoptlib::Problem<T>::TVector;
  using typename cppoptlib::Problem<T>::THessian;

private:
  const std::vector<ChiralityConstraint>& _constraints;
  const DistanceBoundsMatrix& _bounds;

  //! Make an Eigen Vector3d of an atomic index.
  const Eigen::Vector3d _getEigen(const TVector& v, const AtomIndexType& index) {
    assert(v.size() > 3 * index + 2);
    Eigen::Vector3d retv;
    retv << v(3 * index), v(3 * index + 1), v(3 * index + 2);
    return retv;
  }

  T _square(const T& value) { 
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

  const Eigen::Matrix<double, 3, 3> _J(
    const TVector& v,
    const AtomIndexType& i,
    const AtomIndexType& j
  ) {
    Eigen::Matrix<double, 3, 3> retm;

    double a = v(3 * i) - v(3 * j),
           b = v(3 * i + 1) - v(3 * j + 1),
           c = v(3 * i + 2) - v(3 * j + 2);

    retm <<  0, -c,  b,
             c,  0, -a,
            -b,  a,  0;

    return retm;
  }

public:

  DGRefinementProblem(
    const std::vector<ChiralityConstraint>& constraints,
    const DistanceBoundsMatrix& bounds
  ) : 
    _constraints(constraints), 
    _bounds(bounds) {}

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
          distance / (
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
    T sum, intermediate;
    AtomIndexType a, b, c, d;
    T targetVal;

    for(const auto& chiralityConstraint: _constraints) {
      std::tie(a, b, c, d, targetVal) = chiralityConstraint;
      intermediate = targetVal - _getTetrahedronReducedVolume(v, a, b, c, d);
      sum += intermediate * intermediate;
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
      for(AtomIndexType i = 0; i < N; i++) {
        /* Skipping i == alpha is exceptionally important, for the reason that 
         * if i == alpha, although mathematically it would be assumed that
         * _A(i, i) == 0 and _B(i, i) == 0, nevertheless both the upper and 
         * lower bound for this situation are 0, leading to zeroes on 
         * denominators, leading to NaNs!
         */
        if(i == alpha) continue; 
        localGradient += _A(v, i, alpha) + _B(v, alpha, i);
      }

      /* Maybe refactor above into two separate loops, one from 0 -> alpha - 1,
       * then another from alpha + 1 -> N - 1?
       */

      // C 
      for(const auto& chiralityConstraint: _constraints) {
        if(std::get<0>(chiralityConstraint) == alpha) {
          std::tie(std::ignore, a, b, c, target) = chiralityConstraint;

          localGradient += _C(v, alpha, a, b, c, target) * (
            _getEigen(v, a) - _getEigen(v, c)
          ).cross(
            _getEigen(v, b) - _getEigen(v, c)
          );
        } else if(std::get<1>(chiralityConstraint) == alpha) {
          std::tie(a, std::ignore, b, c, target) = chiralityConstraint;

          localGradient += _C(v, a, alpha, b, c, target) * (
            _getEigen(v, a) - _getEigen(v, c)
          ).transpose() * _J(v, b, c);
        } else if(std::get<2>(chiralityConstraint) == alpha) {
          std::tie(a, b, std::ignore, c, target) = chiralityConstraint;

          localGradient += _C(v, a, b, alpha, c, target) * (
            _getEigen(v, a) - _getEigen(v, c)
          ).transpose() * _J(v, b, c);
        } else if(std::get<3>(chiralityConstraint) == alpha) {
          std::tie(a, b, c, std::ignore, target) = chiralityConstraint;

          localGradient += _C(v, a, b, c, alpha, target) * (
            (
              _getEigen(v, a) - _getEigen(v, alpha)
            ).transpose() * (
              _J(v, c, alpha)
              - _J(v, b, alpha) 
            ) - (
              (
                _getEigen(v, b) - _getEigen(v, alpha)
              ).cross(
                _getEigen(v, c) - _getEigen(v, alpha)
              ).transpose()
            )
          );
        }
      }

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
