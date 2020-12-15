/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Class storing atom-pairwise distance bounds
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_DISTANCE_BOUNDS_MATRIX_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_DISTANCE_BOUNDS_MATRIX_H

#include <Eigen/Core>

#include "Molassembler/DistanceGeometry/DistanceGeometry.h"
#include "Molassembler/Modeling/AtomInfo.h"
#include "Molassembler/Conformers.h"

namespace Scine {
namespace Molassembler {

namespace Random {
class Engine;
} // namespace Random

namespace outcome = OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

class DistanceBoundsMatrix {
public:
  static constexpr double defaultLower = 0.0;
  static constexpr double defaultUpper = 100.0;

//!@name Static member functions
//!@{
  static inline double& lowerBound(Eigen::Ref<Eigen::MatrixXd> matrix, const AtomIndex i, const AtomIndex j) {
    if(i < j) {
      return matrix(j, i);
    }

    return matrix(i, j);
  }

  static inline double& upperBound(Eigen::Ref<Eigen::MatrixXd> matrix, const AtomIndex i, const AtomIndex j) {
    if(i < j) {
      return matrix(i, j);
    }

    return matrix(j, i);
  }

  /*! @brief Uses Floyd's algorithm to smooth the matrix
   *
   * @complexity{@math{\Theta(N^3)}}
   */
  static void smooth(Eigen::Ref<Eigen::MatrixXd> matrix);
//!@}

//!@name Member types
//!@{
  using BoundsMatrix = Eigen::MatrixXd;
//!@}

//!@name Special member functions
//!@{
  //! Initializes nothing
  DistanceBoundsMatrix();

  //! Initializes the matrix with defaultLower and defaultUpper
  explicit DistanceBoundsMatrix(long unsigned N);

  //! Sets the underlying matrix as the argument
  explicit DistanceBoundsMatrix(Eigen::MatrixXd matrix);

  /**
   * @brief Constructs a bounds matrix from a graph and a bounds matrix
   *
   * Replaces default lower bounds with the sum of the elements' vdw radii.
   */
  template<typename PrivateGraph>
  DistanceBoundsMatrix(
    const PrivateGraph& inner,
    BoundsMatrix bounds
  ) : matrix_(std::move(bounds)) {
    // Populate the lower bounds if no explicit information is present
    const AtomIndex N = inner.V();
    for(AtomIndex i = 0; i < N - 1; ++i) {
      for(AtomIndex j = i + 1; j < N; ++j) {
        // i < j in all cases -> lower bound is (j, i), upper bound is (i, j)
        if(matrix_(j, i) == 0) {
          matrix_(j, i) = (
            AtomInfo::vdwRadius(inner.elementType(i))
            + AtomInfo::vdwRadius(inner.elementType(j))
          );
        }

        // Ensure upper bound has a default value for Floyd's algorithm
        if(matrix_(i, j) == 0.0) {
          matrix_(i, j) = defaultUpper;
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

  //! Smoothes the underlying matrix using smooth(Eigen::Ref<Eigen::MatrixXd>)
  void smooth();
//!@}

//!@name Information
//!@{
  /*! @brief Access to upper bound of unordered indices
   *
   * @complexity{@math{\Theta(1)}}
   */
  inline double upperBound(AtomIndex i, AtomIndex j) const {
    if(i < j) {
      return matrix_(i, j);
    }

    return matrix_(j, i);
  }

  /*! @brief Access to lower bound of unordered indices
   *
   * @complexity{@math{\Theta(1)}}
   */
  inline double lowerBound(AtomIndex i, AtomIndex j) const {
    if(i < j) {
      return matrix_(j, i);
    }

    return matrix_(i, j);
  }

  /*! @brief Checks for cases in which the lower bound is greater than the upper bound
   *
   * @complexity{@math{\Theta(N^2)}}
   */
  unsigned boundInconsistencies() const;

  /** @brief Nonmodifiable access to bounds matrix
   *
   * @complexity{@math{\Theta(1)}}
   */
  const Eigen::MatrixXd& access() const;

  /*! @brief Generate a distance matrix
   *
   * Generates a distance matrix from the contained bounds without altering
   * underlying state.
   *
   * @complexity{@math{O(N^5)}}
   */
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(
    Random::Engine& engine,
    Partiality partiality = Partiality::All
  ) const noexcept;

  /*! @brief Squares all bounds and returns a new matrix
   *
   * @complexity{@math{\Theta(N^2)}}
   */
  Eigen::MatrixXd makeSquaredBoundsMatrix() const;

  /*! @brief Yields the number of particles
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned N() const;
//!@}

private:
  Eigen::MatrixXd matrix_;

  inline double& lowerBound_(const AtomIndex i, const AtomIndex j) {
    return matrix_(
      std::max(i, j),
      std::min(i, j)
    );
  }

  inline double& upperBound_(const AtomIndex i, const AtomIndex j) {
    return matrix_(
      std::min(i, j),
      std::max(i, j)
    );
  }

};

} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine

#endif
