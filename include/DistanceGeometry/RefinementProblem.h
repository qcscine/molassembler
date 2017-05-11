#ifndef INCLUDE_DG_REFINEMENT_PROBLEM_H
#define INCLUDE_DG_REFINEMENT_PROBLEM_H

#include "cppoptlib/problem.h"
#include "common_typedefs.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"

#include <Eigen/Dense>

/* NOTES
 * - RefinementProblem has two full gradient implementations, one which is 
 *   properly separated into subproblems and is written in a very legible
 *   style. The other is optimized as best as possible and has hence suffered
 *   some legibility impairments. In testing, the two are compared to ensure
 *   that no errors are introduced by any further optimizations. 
 * - Profiling some DG runs leads to the conclusion that most time is spent in
 *   value() instead of gradient() so putting some effort into that is probably 
 *   most valuable.
 */

/* TODO
 * - Optimize 
 */

namespace MoleculeManip {

namespace DistanceGeometry {

class RefinementProblem : public cppoptlib::Problem<double> {
public:
/* Typedefs */
  using typename cppoptlib::Problem<double>::TVector;
  using typename cppoptlib::Problem<double>::THessian;

private:
/* Private member functions */
  //! Make an Eigen Vector4d of an atomic index.
  template<unsigned vectorSize = 4>
  inline auto _getPos(const TVector& v, const AtomIndexType& index) const {
    static_assert(
      vectorSize == 3 || vectorSize == 4,
      "vectorSize in _getPos instantiation must be 3 or 4"
    );
    assert(v.size() > static_cast<long>(4 * index + 3));

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

  inline double _square(const double& value) const { 
    return value * value;
  }

  /*!
   * Reduced since the factor 6 is not multiplied here to the result of the
   * vector dot product.
   */
  inline double _getTetrahedronReducedVolume(
    const TVector& v,
    const std::array<AtomIndexType, 4>& indices
  ) const {
    const Eigen::Vector3d vecL = _getPos<3>(v, indices.at(3));

    return (
      _getPos<3>(v, indices.at(0)) - vecL
    ).dot(
      (
        _getPos<3>(v, indices.at(1)) - vecL
      ).cross(
        _getPos<3>(v, indices.at(2)) - vecL
      )
    );
  }


public:
/* Typedefs */
  using CallbackFunction = std::function<void(const cppoptlib::Criteria<double>&, const TVector&, const RefinementProblem&)>;

/* State */
  // More or less constant
  const std::vector<ChiralityConstraint>& constraints;
  const Eigen::MatrixXd squaredBounds;
  CallbackFunction callbackClosure;
  const unsigned N;

  // True mutable state
  mutable bool compress = false;
  mutable double fourthDimWeight = 0.0;

/* Constructors */
  RefinementProblem(
    const std::vector<ChiralityConstraint>& constraints,
    const DistanceBoundsMatrix& bounds
  );

  RefinementProblem(
    const std::vector<ChiralityConstraint>& constraints,
    const DistanceBoundsMatrix& bounds,
    const CallbackFunction& callback
  );

/* Information */
  bool callback(
    const cppoptlib::Criteria<double>& state,
    const TVector& x
  ) const;

  double chiralError(const TVector& v) const;

  double extraDimensionError(const TVector& v) const;

  double distanceError(const TVector& v) const;

  /*!
   * Required for cppoptlib: Calculates gradient for specified coordinates.
   * Optimized version of gradient implementation.
   */
  void gradient(const TVector& v, TVector& grad) const;

  /*! Properly separated and legible reference gradient implementation.
   * Due to intensive testing, this is presumed correct. If you want to speed
   * up DG, put your work into the optimized version.
   */
  void referenceGradient(const TVector& v, TVector& gradient) const;

  void referenceGradientA(const TVector& v, TVector& gradient) const;
  void referenceGradientB(const TVector& v, TVector& gradient) const;
  void referenceGradientC(const TVector& v, TVector& gradient) const;
  void referenceGradientD(const TVector& v, TVector& gradient) const;

  // Helper functions for reference implementation readability
  inline const double& upperBound(const AtomIndexType& i, const AtomIndexType& j) const {
    return squaredBounds(
      std::min(i, j),
      std::max(i, j)
    );
  }

  inline const double& lowerBound(const AtomIndexType& i, const AtomIndexType& j) const {
    return squaredBounds(
      std::max(i, j),
      std::min(i, j)
    );
  }

  //! Inverts all y coordinates of the passed positions vector
  void invertY(TVector& v) const;

  /*! Calculates the fraction of correct chirality constraints over the total
   * amount of non-zero chirality constraints. (Planar chirality constraints are
   * not regarded as they cannot be the wrong "sign".)
   */
  double proportionCorrectChiralityConstraints(const TVector& v) const;

  //! Calculates error function value for specified coordinates.
  double value(const TVector& v) override;
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
