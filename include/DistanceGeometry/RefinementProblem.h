#ifndef INCLUDE_DG_REFINEMENT_PROBLEM_H
#define INCLUDE_DG_REFINEMENT_PROBLEM_H

#include "cppoptlib/problem.h"
#include "common_typedefs.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"

#include <Eigen/Dense>

/* TODO
 * - Optimize 
 *   - ... ?
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
  const std::vector<ChiralityConstraint>& constraints;
  const BoundsMatrix squaredBounds;
  CallbackFunction callbackClosure;
  mutable bool compress = false;

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

  void alternateGradient(const TVector& v, TVector& gradient) const;

  /*!
   * Required for cppoptlib: Calculates gradient for specified coordinates.
   */
  void gradient(const TVector& v, TVector& grad) const;

  void gradientA(const TVector& v, TVector& gradient) const;

  void gradientB(const TVector& v, TVector& gradient) const;

  void gradientC(const TVector& v, TVector& gradient) const;

  void gradientD(const TVector& v, TVector& gradient) const;

  //!  First term of gradient expansion, used within gradient()
  inline Eigen::Vector4d gradientTermA(
    const TVector& v,
    const AtomIndexType& i, 
    const AtomIndexType& j
  ) const {
    const Eigen::Vector4d retv = _getPos(v, j) - _getPos(v, i);

    return retv * 4 * (
      retv.squaredNorm() / squaredBounds.upperBound(i, j)
      - 1
    ) / squaredBounds.upperBound(i, j);
  }

  //! Second term of gradient expansion, used within gradient()
  inline Eigen::Vector4d gradientTermB(
    const TVector& v,
    const AtomIndexType& i,
    const AtomIndexType& j
  ) const {
    const Eigen::Vector4d retv = _getPos(v, j) - _getPos(v, i);

    return retv * 8 * squaredBounds.lowerBound(i, j) * (
      squaredBounds.lowerBound(i, j) - retv.squaredNorm()
    ) / std::pow(
      squaredBounds.lowerBound(i, j) + retv.squaredNorm(),
      3
    );
  }

  //! Third term of gradient expansion, used within gradient()
  inline double gradientTermC(
    const TVector& v,
    const std::array<AtomIndexType, 4>& indices,
    const double& target
  ) const {
    return (-2) * (
      target - _getTetrahedronReducedVolume(v, indices)
    );
  }

  void invertY(TVector& v) const;

  double proportionCorrectChiralityConstraints(const TVector& v) const;

  /*! 
   * Calculates value for specified coordinates.
   */
  double value(const TVector& v) override;
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
