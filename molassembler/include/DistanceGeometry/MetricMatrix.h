#ifndef INCLUDE_DG_METRIC_MATRIX_HPP
#define INCLUDE_DG_METRIC_MATRIX_HPP

#include <Eigen/Core>

#include "DistanceGeometry/DistanceGeometry.h"

/*! @file
 *
 * In the Distance Geometry algorithm, a metric matrix is generated from the
 * distance bounds. This file hold the declaration of the implementation of that
 * matrix.
 */

namespace molassembler {

namespace DistanceGeometry {

class MetricMatrix {
private:
/* Underlying matrix representation */
  Eigen::MatrixXd _matrix;

  void _constructFromTemporary(Eigen::MatrixXd&& distances);

public:
/* Constructors */
  MetricMatrix() = delete;
  explicit MetricMatrix(Eigen::MatrixXd&& matrix); 
  explicit MetricMatrix(const Eigen::MatrixXd& matrix); 

/* Information */
  //! Allow const ref access to underlying matrix
  const Eigen::MatrixXd& access() const;

  /*!
   * Embeds itself into 4D space, returning a dynamically sized Matrix where
   * every column vector is the coordinates of a particle
   */
  Eigen::MatrixXd embed() const;

/* Operators */
  bool operator == (const MetricMatrix& other) const;
  friend std::ostream& operator << (std::ostream& os, const MetricMatrix& metricMatrix);
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
