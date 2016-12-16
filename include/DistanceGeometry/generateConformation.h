#ifndef INCLUDE_DG_GENERATE_CONFORMATION_H
#define INCLUDE_DG_GENERATE_CONFORMATION_H

#include "Molecule.h"
#include "Types/PositionCollection.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/MetricMatrix.h"
#include "DistanceGeometry/DGRefinementProblem.h"

#include "cppoptlib/meta.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"

#include <vector>
#include <Eigen/Core>

/* TODO
 * - consider a refactor into a class with various member functions instead of
 *   so many "free-roaming" functions. Would alleviate a lot of object passing.
 */

namespace MoleculeManip {

namespace DistanceGeometry {

Eigen::Vector3d getPos(
  const Eigen::MatrixXd& positions,
  const AtomIndexType& index
) {
  Eigen::Vector3d retv;
  retv = positions.col(index);
  return retv;
}

double evaluateChiralityConstraint(
  const ChiralityConstraint& chiralityConstraint,
  const Eigen::MatrixXd& positions
) {
  AtomIndexType i, j, k, l;
  std::tie(i, j, k, l, std::ignore) = chiralityConstraint;

  // V = (1 - 4) * [ (2 - 4) x (3 - 4) ]
  return (
    (
      getPos(positions, i)
      - getPos(positions, l)
    ).dot(
      (
       getPos(positions, j)
       - getPos(positions, l)
      ).cross(
        getPos(positions, k)
        - getPos(positions, l)
      )
    )
  );
}

bool moreThanHalfChiralityConstraintsIncorrect(
  const Eigen::MatrixXd& positions,
  const std::vector<ChiralityConstraint>& chiralityConstraints
) {
  unsigned totalNonZeroConstraints = 0, incorrectNonZeroConstraints = 0;
  for(const auto& chiralityConstraint : chiralityConstraints) {
    auto& target = std::get<4>(chiralityConstraint);

    if(target != 0.0) {
      totalNonZeroConstraints += 1;
    }

    auto eval = evaluateChiralityConstraint(
      chiralityConstraint,
      positions
    );

    if( // can this be simplified? -> sign bit XOR?
      ( eval < 0 && target > 0)
      || (eval > 0 && target < 0)
    ) {
      incorrectNonZeroConstraints += 1;
    }
  }

  return (
    incorrectNonZeroConstraints / totalNonZeroConstraints
    > 0.5
  );

}

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
) {
  auto distanceBoundsMatrix = molecule.getDistanceBoundsMatrix();
  auto chiralityConstraints = molecule.getChiralityConstraints();

  Eigen::MatrixXd embeddedPositions;
  bool acceptableConformation;

  DGRefinementProblem<double> problem(
    chiralityConstraints,
    distanceBoundsMatrix
  );

  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGConjugatedGradientDescentSolver;

  do {

    MetricMatrix metric(
      distanceBoundsMatrix.generateDistanceMatrix(
        metrization
      )
    );

    embeddedPositions = metric.embed(embedding);

    if(moreThanHalfChiralityConstraintsIncorrect(
      embeddedPositions,
      chiralityConstraints
    )) {
      embeddedPositions.row(2) *= -1;
    }

    Eigen::VectorXd vectorizedPositions(
      Eigen::Map<Eigen::VectorXd>(
        embeddedPositions.data(),
        embeddedPositions.cols() * embeddedPositions.rows()
      )
    );

    DGConjugatedGradientDescentSolver.minimize(problem, vectorizedPositions);

    // TODO TEST acceptable or not -> acceptableConformation

  } while(!acceptableConformation);

  /* somehow convert Eigen matrix into PositionCollection */
  Delib::PositionCollection positions;
  // for every column in embedded?
  //   positions.push_back(Delib::Position([Eigen::Vector3d] column))
  return positions;
}



} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
