#include <vector>
#include <Eigen/Core>

#include "Types/PositionCollection.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"

namespace MoleculeManip {

namespace DistanceGeometry {

bool refine(
  Eigen::MatrixXd& embedded,
  const EmbeddingOption& embedding,
  const MoleculeManip::Molecule& molecule
);

Delib::PositionCollection generateConformation(
  const MoleculeManip::Molecule& molecule,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
) {
  DistanceBoundsMatrix distanceBounds(molecule);
  Eigen::MatrixXd embedded;
  bool acceptableConformation;

  do {
    MetricMatrix metric = distanceBounds.toMetricMatrix(metrization);
    embedded = metric.embed(embedding);

    acceptableConformation = refine(embedded, embedding, molecule);
  } while(!acceptableConformation);

  /* somehow convert Eigen matrix into PositionCollection */
  Delib::PositionCollection positions;
  // for every column in embedded?
  //   positions.push_back(Delib::Position([Eigen::Vector3d] column))
  return positions;
}

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip
