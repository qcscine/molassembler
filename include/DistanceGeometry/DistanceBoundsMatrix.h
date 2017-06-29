#ifndef INCLUDE_DG_DISTANCE_BOUNDS_MATRIX_HPP
#define INCLUDE_DG_DISTANCE_BOUNDS_MATRIX_HPP

#include "DistanceGeometry/DistanceGeometry.h"
#include "common_typedefs.h"

#include <Eigen/Core>
#include <random>

namespace MoleculeManip {

namespace DistanceGeometry {

struct BoundsMatrix {
/* State */
  Eigen::MatrixXd matrix;
  
/* Constructors */
  explicit BoundsMatrix(const unsigned& N) {
    matrix.resize(N, N);
    matrix.setZero();
  }

  explicit BoundsMatrix(Eigen::MatrixXd passMatrix) : matrix(passMatrix) {}

/* Modification */
  double& lowerBound(const unsigned& i, const unsigned& j) {
    return matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  double& upperBound(const unsigned& i, const unsigned& j) {
    return matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }

/* Smoothing algorithms */
  //! Smooth the bounds once with the triangle inequality
  bool triangleInequalitySmooth() {
    /* Floyd's algorithm
     * Could be refactored slightly that when something is changed, these loops
     * exit and identical loops without the bool setter are run (minimal speed gains)
     */
    bool changedSomething = false;
    const unsigned N = matrix.cols();

    for(unsigned i = 0; i < N - 1; i++) {
      for(unsigned j = i + 1; j < N; j++) {
        for(unsigned k = 0; k < N; k++) {
          if(k == i || k == j) {
            continue;
          }

          if(upperBound(i, j) > upperBound(i, k) + upperBound(k, j)) {
            upperBound(i, j) = upperBound(i, k) + upperBound(k, j);
            changedSomething = true;
          }

          if(lowerBound(i, j) < lowerBound(i, k) - upperBound(k, j)) {
            lowerBound(i, j) = lowerBound(i, k) - upperBound(k, j);
            changedSomething = true;
          } else if(lowerBound(i, j) < lowerBound(j, k) - upperBound(k, i)) {
            lowerBound(i, j) = lowerBound(j, k) - upperBound(k, i);
            changedSomething = true;
          }

          assert(lowerBound(i, j) <= upperBound(i, j));
        }
      }
    }

    return changedSomething;
  }

  void smooth() {
    const unsigned maxIter = 100;
    for(
      unsigned iter = 0;
      ( // run as long as
        iter < maxIter // we do not exceed 100 iterations
        && triangleInequalitySmooth() // and smoothing changes something
      ); 
      iter++
    ) {
      continue;
    }
  }

/* Information */
  double lowerBound(const unsigned& i, const unsigned& j) const {
    return matrix(
      std::max(i, j),
      std::min(i, j)
    );
  }

  double upperBound(const unsigned& i, const unsigned& j) const {
    return matrix(
      std::min(i, j),
      std::max(i, j)
    );
  }
};

class DistanceBoundsMatrix {
private:
/* Data members */
  // Matrix data
  BoundsMatrix _boundsMatrix;

  // Randomness state
  mutable std::mt19937 _randomEngine;

/* Private members */
  //! Constructor helper, initializes randomness state
  void _initRandomEngine();

public:
  const unsigned N;

/* Constructors */
  DistanceBoundsMatrix() = delete;
  explicit DistanceBoundsMatrix(const unsigned& N) noexcept;
  explicit DistanceBoundsMatrix(const Eigen::MatrixXd& matrix) noexcept;

  DistanceBoundsMatrix(const DistanceBoundsMatrix& other) noexcept;
  DistanceBoundsMatrix(DistanceBoundsMatrix&& other) noexcept;
  DistanceBoundsMatrix& operator = (const DistanceBoundsMatrix& other) noexcept;
  DistanceBoundsMatrix& operator = (DistanceBoundsMatrix&& other) noexcept;

/* Modifiers */
  //! Smooth until the matrix does not change
  void smooth();

  bool setLowerBound(
    const unsigned& i,
    const unsigned& j,
    const double& newLowerBound
  );

  bool setUpperBound(
    const unsigned& i,
    const unsigned& j,
    const double& newUpperBound
  );

/* Information */
  //! Access the underlying matrix
  const Eigen::MatrixXd& access() const;

  /*! 
   * Returns a count of instances where lower and upper bounds are logically 
   * wrong, i.e. when the lower bound for a pair of atoms is larger than the 
   * upper bound.
   */
  unsigned boundInconsistencies() const;

  /*!
   * Returns a distance matrix with randomly chosen distances between the
   * bounds
   */
  Eigen::MatrixXd generateDistanceMatrix(
    const MetrizationOption& metrization = MetrizationOption::full
  ) const;

  //! Access a lower bound
  double lowerBound(
    const unsigned& i,
    const unsigned& j
  ) const;

  //! Returns a bounds matrix with the current bounds squared
  Eigen::MatrixXd makeSquaredBoundsMatrix() const;

  //! Access an upper bound
  double upperBound(
    const unsigned& i,
    const unsigned& j
  ) const;
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
