#ifndef INCLUDE_DG_DISTANCE_BOUNDS_MATRIX_HPP
#define INCLUDE_DG_DISTANCE_BOUNDS_MATRIX_HPP

#include <Eigen/Core>

#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/MetricMatrix.h"

namespace MoleculeManip {

namespace DistanceGeometry {

class DistanceBoundsMatrix {
private:
  /* Underlying matrix representation */
  Eigen::MatrixXd _matrix;

  void _setUpperBound(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const double& value
  );

  void _setLowerBound(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const double& value
  );

public:
  /* Constructors */
  DistanceBoundsMatrix() = delete;
  DistanceBoundsMatrix(
    const MoleculeManip::Molecule& molecule
  );

  MetricMatrix toMetricMatrix(
    const MetrizationOption& metrization
  );

};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
