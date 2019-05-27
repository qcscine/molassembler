/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Class storing atom-pairwise distance bounds
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_DISTANCE_BOUNDS_MATRIX_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_DISTANCE_BOUNDS_MATRIX_H

#include "boost_outcome/outcome.hpp"
#include <Eigen/Core>
#include "temple/Random.h"

#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/Modeling/AtomInfo.h"
#include "molassembler/Conformers.h"

namespace Scine {

namespace molassembler {

namespace random {
class Engine;
} // namespace random

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

class DistanceBoundsMatrix {
private:
  inline double& _lowerBound(const AtomIndex i, const AtomIndex j) {
    return _matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  inline double& _upperBound(const AtomIndex i, const AtomIndex j) {
    return _matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

public:
  static constexpr double defaultLower = 0.0;
  static constexpr double defaultUpper = 100.0;

//!@name Static member functions
//!@{
  static inline double& lowerBound(Eigen::MatrixXd& matrix, const AtomIndex i, const AtomIndex j) {
    return matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  static inline double& upperBound(Eigen::MatrixXd& matrix, const AtomIndex i, const AtomIndex j) {
    return matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  static inline double lowerBound(const Eigen::MatrixXd& matrix, const AtomIndex i, const AtomIndex j) {
    return matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  static inline double upperBound(const Eigen::MatrixXd& matrix, const AtomIndex i, const AtomIndex j) {
    return matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  static void smooth(Eigen::MatrixXd& matrix);
//!@}

//!@name Member types
//!@{
  using BoundsMatrix = Eigen::MatrixXd;
//!@}

//!@name Special member functions
//!@{
  DistanceBoundsMatrix();

  explicit DistanceBoundsMatrix(long unsigned N);

  explicit DistanceBoundsMatrix(Eigen::MatrixXd matrix);

  template<typename InnerGraph>
  DistanceBoundsMatrix(
    const InnerGraph& inner,
    BoundsMatrix bounds
  ) : _matrix(std::move(bounds)) {
    // Populate the lower bounds if no explicit information is present
    const AtomIndex N = inner.N();
    for(AtomIndex i = 0; i < N - 1; ++i) {
      for(AtomIndex j = i + 1; j < N; ++j) {
        // i < j in all cases -> lower bound is (j, i), upper bound is (i, j)
        if(_matrix(j, i) == 0) {
          _matrix(j, i) = (
            AtomInfo::vdwRadius(inner.elementType(i))
            + AtomInfo::vdwRadius(inner.elementType(j))
          );
        }

        // Ensure upper bound has a default value for Floyd's algorithm
        if(_matrix(i, j) == 0.0) {
          _matrix(i, j) = defaultUpper;
        }
      }
    }

    assert(boundInconsistencies() == 0);
  }
//!@}

//!@name Modifiers
//!@{
  bool setUpperBound(AtomIndex i, AtomIndex j, double newUpperBound);
  bool setLowerBound(AtomIndex i, AtomIndex j, double newLowerBound);

  void smooth();
//!@}

//!@name Information
//!@{
  inline double upperBound(AtomIndex i, AtomIndex j) const {
    return _matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  inline double lowerBound(AtomIndex i, AtomIndex j) const {
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
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(random::Engine& engine, Partiality partiality = Partiality::All) const noexcept;

  Eigen::MatrixXd makeSquaredBoundsMatrix() const;

  unsigned N() const;
//!@}

private:
  Eigen::MatrixXd _matrix;
};

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine

#endif
