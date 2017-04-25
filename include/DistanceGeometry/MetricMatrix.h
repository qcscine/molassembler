#ifndef INCLUDE_DG_METRIC_MATRIX_HPP
#define INCLUDE_DG_METRIC_MATRIX_HPP

#include <Eigen/Core>

#include "DistanceGeometry/DistanceGeometry.h"

namespace MoleculeManip {

namespace DistanceGeometry {

class MetricMatrix {
private:
/* Underlying matrix representation */
  Eigen::MatrixXd _matrix;

  void _constructFromTemporary(Eigen::MatrixXd&& matrix);

public:
/* Constructors */
  MetricMatrix() = delete;
  MetricMatrix(Eigen::MatrixXd&& matrix); // want to be able to modify temporary
  MetricMatrix(const Eigen::MatrixXd& matrix); 

/* Information */
  //! Allow const ref access to underlying matrix
  const Eigen::MatrixXd& access() const;

  /*!
   * Embeds itself into 3D or 4D space depending on the embedding option,
   * returning a dynamically sized Matrix where every column vector is the
   * coordinates of a particle
   */
  Eigen::MatrixXd embed(const EmbeddingOption& embedding) const;

/* Operators */
  bool operator == (const MetricMatrix& other) const;
  friend std::ostream& operator << (std::ostream& os, const MetricMatrix& matrix);
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
