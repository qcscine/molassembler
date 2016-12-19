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
  mutable std::mt19937 _randomEngine;

  /* Constructor helper */
  void _initRandomEngine();

public:
  enum class SmoothingAlgorithm {
    Naive,
    Custom
  };

  /* Constructors */
  DistanceBoundsMatrix() = delete;
  DistanceBoundsMatrix(const unsigned& N);
  DistanceBoundsMatrix(Eigen::MatrixXd matrix);

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

  void triangleInequalitySmooth(const SmoothingAlgorithm& algorithmChoice);

  /*!
   * Returns a distance matrix with randomly chosen distances between the
   * bounds
   */
  Eigen::MatrixXd generateDistanceMatrix(
    const MetrizationOption& metrization
  ) const;
};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
