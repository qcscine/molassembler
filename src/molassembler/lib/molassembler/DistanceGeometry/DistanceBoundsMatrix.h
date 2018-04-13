#ifndef INCLUDE_DG_DISTANCE_BOUNDS_MATRIX_HPP
#define INCLUDE_DG_DISTANCE_BOUNDS_MATRIX_HPP

#include "boost_outcome/outcome.hpp"
#include "DistanceGeometry/DistanceGeometry.h"
#include "AtomInfo.h"

#include <Eigen/Core>
#include "temple/Random.h"

/*! @file
 *
 * Contains the implementation of a class that stores distance bounds.
 */

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

class DistanceBoundsMatrix {
private:
  Eigen::MatrixXd _matrix;

  inline double& _lowerBound(const AtomIndexType& i, const AtomIndexType& j) {
    return _matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  inline double& _upperBound(const AtomIndexType& i, const AtomIndexType& j) {
    return _matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

public:
  static inline double& lowerBound(Eigen::MatrixXd& matrix, const AtomIndexType& i, const AtomIndexType& j) {
    return matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  static inline double& upperBound(Eigen::MatrixXd& matrix, const AtomIndexType& i, const AtomIndexType& j) {
    return matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  static inline double lowerBound(const Eigen::MatrixXd& matrix, const AtomIndexType& i, const AtomIndexType& j) {
    return matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  static inline double upperBound(const Eigen::MatrixXd& matrix, const AtomIndexType& i, const AtomIndexType& j) {
    return matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  using BoundList = std::vector<
    std::tuple<AtomIndexType, AtomIndexType, ValueBounds>
  >;

  DistanceBoundsMatrix();

  explicit DistanceBoundsMatrix(const unsigned& N);

  explicit DistanceBoundsMatrix(const Eigen::MatrixXd& matrix);

  template<typename Molecule>
  DistanceBoundsMatrix(const Molecule& molecule, const BoundList& bounds) {
    unsigned N = molecule.numAtoms();
    _matrix.resize(N, N);
    _matrix.triangularView<Eigen::Lower>().setZero();
    _matrix.triangularView<Eigen::StrictlyUpper>().setConstant(100);

    // Populate the matrix with explicit bounds
    AtomIndexType a, b;
    ValueBounds bound;
    for(const auto& boundTuple : bounds) {
      std::tie(a, b, bound) = boundTuple;

      setUpperBound(a, b, bound.upper);
      setLowerBound(a, b, bound.lower);
    }

    // Populate the lower bounds if no explicit information is present
    for(AtomIndexType i = 0; i < N - 1; ++i) {
      for(AtomIndexType j = i + 1; j < N; ++j) {
        /* setting the bounds will fail for bonded pairs as those have strict
         * bounds already and the fairly high sum of vdw would lead to
         * inconsistencies
         */
        if(lowerBound(i, j) == 0) {
          setLowerBound(
            i,
            j,
            AtomInfo::vdwRadius(
              molecule.getElementType(i)
            ) + AtomInfo::vdwRadius(
              molecule.getElementType(j)
            )
          );
        }
      }
    }

    assert(boundInconsistencies() == 0);
  }

  bool setUpperBound(const AtomIndexType& i, const AtomIndexType& j, double newUpperBound);

  bool setLowerBound(const AtomIndexType& i, const AtomIndexType& j, const double& newLowerBound);

  static void smooth(Eigen::MatrixXd& matrix);

  void smooth();

  inline double upperBound(const AtomIndexType& i, const AtomIndexType& j) const {
    return _matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  inline double lowerBound(const AtomIndexType& i, const AtomIndexType& j) const {
    return _matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  unsigned boundInconsistencies() const;

  const Eigen::MatrixXd& access() const;

  /* Generates a distance matrix from the contained bounds without destroying
   * the inherent state.
   *
   * Allocates another N*N double matrix. When resource constrained, this is
   * not a good idea.
   */
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(Partiality partiality = Partiality::All) const noexcept;

  Eigen::MatrixXd makeSquaredBoundsMatrix() const;

  unsigned N() const;
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
