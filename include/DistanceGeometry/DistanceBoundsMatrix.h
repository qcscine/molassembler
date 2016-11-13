#ifndef INCLUDE_DG_DISTANCE_BOUNDS_MATRIX_HPP
#define INCLUDE_DG_DISTANCE_BOUNDS_MATRIX_HPP

#include "DistanceGeometry/DistanceGeometry.h"
#include "common_typedefs.h"

#include <Eigen/Core>
#include <random>

namespace MoleculeManip {

namespace DistanceGeometry {

class DistanceBoundsMatrix {
private:
  /* Underlying matrix representation */
  Eigen::MatrixXd _matrix;
  const unsigned _N;
  std::vector<unsigned> _seeds;
  std::seed_seq _seedSequence;
  std::mt19937 _randomEngine;

public:
  /* Constructors */
  DistanceBoundsMatrix() = delete;
  DistanceBoundsMatrix(const unsigned& N);

  /* Accessor & Modifier */
  decltype(_matrix(1, 2))& upperBound(
    const unsigned& i,
    const unsigned& j
  );
  decltype(_matrix(1, 2))& lowerBound(
    const unsigned& i,
    const unsigned& j
  );
  double upperBound(
    const unsigned& i,
    const unsigned& j
  ) const;
  double lowerBound(
    const unsigned& i,
    const unsigned& j
  ) const;

  void processDistanceConstraints(
    const std::vector<DistanceConstraint>& constraints
  );

  /*!
   * Returns a distance matrix with randomly chosen distances between the
   * bounds
   */
  Eigen::MatrixXd generateDistanceMatrix(
    const MetrizationOption& metrization
  );
};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
