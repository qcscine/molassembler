#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_DISTANCE_BOUNDS_MATRIX_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_DISTANCE_BOUNDS_MATRIX_H

#include "boost_outcome/outcome.hpp"
#include <Eigen/Core>
#include "temple/Random.h"

#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/AtomInfo.h"


/*! @file
 *
 * Contains the implementation of a class that stores distance bounds.
 */

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

namespace DistanceGeometry {

class DistanceBoundsMatrix {
private:
  inline double& _lowerBound(const AtomIndexType i, const AtomIndexType j) {
    return _matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  inline double& _upperBound(const AtomIndexType i, const AtomIndexType j) {
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
  static inline double& lowerBound(Eigen::MatrixXd& matrix, const AtomIndexType i, const AtomIndexType j) {
    return matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  static inline double& upperBound(Eigen::MatrixXd& matrix, const AtomIndexType i, const AtomIndexType j) {
    return matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  static inline double lowerBound(const Eigen::MatrixXd& matrix, const AtomIndexType i, const AtomIndexType j) {
    return matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  static inline double upperBound(const Eigen::MatrixXd& matrix, const AtomIndexType i, const AtomIndexType j) {
    return matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  static void smooth(Eigen::MatrixXd& matrix);
//!@}

//!@name Member types
//!@{
  using BoundsList = std::map<
    std::array<AtomIndexType, 2>,
    ValueBounds
  >;
//!@}

//!@name Special member functions
//!@{
  DistanceBoundsMatrix();

  explicit DistanceBoundsMatrix(unsigned N);

  explicit DistanceBoundsMatrix(Eigen::MatrixXd matrix);

  template<typename Molecule>
  DistanceBoundsMatrix(
    const Molecule& molecule,
    const BoundsList& bounds
  ) : DistanceBoundsMatrix {molecule.numAtoms()}
  {
    // Populate the matrix with explicit bounds
    for(const auto& mapPair : bounds) {
      setLowerBound(mapPair.first.front(), mapPair.first.back(), mapPair.second.lower);
      setUpperBound(mapPair.first.front(), mapPair.first.back(), mapPair.second.upper);
    }

    // Populate the lower bounds if no explicit information is present
    const AtomIndexType N = molecule.numAtoms();
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
//!@}

//!@name Modifiers
//!@{
  bool setUpperBound(AtomIndexType i, AtomIndexType j, double newUpperBound);
  bool setLowerBound(AtomIndexType i, AtomIndexType j, double newLowerBound);

  void smooth();
//!@}

//!@name Information
//!@{
  inline double upperBound(AtomIndexType i, AtomIndexType j) const {
    return _matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

  inline double lowerBound(AtomIndexType i, AtomIndexType j) const {
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
//!@}

private:
  Eigen::MatrixXd _matrix;
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
