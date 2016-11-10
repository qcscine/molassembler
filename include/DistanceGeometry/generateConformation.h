#include "Molecule.h"
#include "Types/PositionCollection.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/MetricMatrix.h"

#include <vector>
#include <Eigen/Core>

namespace MoleculeManip {

namespace DistanceGeometry {

bool refine(
  Eigen::MatrixXd& embedded,
  const EmbeddingOption& embedding,
  const Molecule& molecule
);

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
) {
  auto distanceBoundsMatrix = molecule.getDistanceBoundsMatrix();
  Eigen::MatrixXd embedded;
  bool acceptableConformation;

  do {
    MetricMatrix metric(
      distanceBoundsMatrix.generateDistanceMatrix(
        metrization
      )
    );
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
