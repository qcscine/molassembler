#ifndef INCLUDE_DG_METRIC_MATRIX_HPP
#define INCLUDE_DG_METRIC_MATRIX_HPP

#include <Eigen/Core>

#include "DistanceGeometry/DistanceGeometry.h"
#include "Molecule.h"

namespace MoleculeManip {

namespace DistanceGeometry {

class MetricMatrix {
private:
  /* Underlying matrix representation */
  Eigen::MatrixXd matrix;

public:
  /* Constructors */
  MetricMatrix() = delete;
  MetricMatrix(
    const Eigen::MatrixXd& matrix
  );

  /*!
   * Embeds itself into 3D or 4D space depending on the embedding option,
   * returning a dynamically sized Matrix
   */
  Eigen::MatrixXd embed(
    const EmbeddingOption& embedding
  );
};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
