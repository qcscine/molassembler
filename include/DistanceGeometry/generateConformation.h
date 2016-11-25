#include "Molecule.h"
#include "Types/PositionCollection.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/MetricMatrix.h"

#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"

#include <vector>
#include <Eigen/Core>

/* TODO
 * - consider a refactor into a class with various member functions instead of
 *   so many "free-roaming" functions. Would alleviate a lot of object passing.
 */

namespace MoleculeManip {

namespace DistanceGeometry {

bool refine(
  Eigen::MatrixXd& embedded,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  DistanceBoundsMatrix& distanceBoundsMatrix
);

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
) {
  auto distanceBoundsMatrix = molecule.getDistanceBoundsMatrix();
  Eigen::MatrixXd embeddedPositions;
  bool acceptableConformation;

  do {
    MetricMatrix metric(
      distanceBoundsMatrix.generateDistanceMatrix(
        metrization
      )
    );

    embeddedPositions = metric.embed(embedding);

    acceptableConformation = refine(
      embeddedPositions,
      molecule.getChiralityConstraints(),
      distanceBoundsMatrix
    );

  } while(!acceptableConformation);

  /* somehow convert Eigen matrix into PositionCollection */
  Delib::PositionCollection positions;
  // for every column in embedded?
  //   positions.push_back(Delib::Position([Eigen::Vector3d] column))
  return positions;
}

double evaluateChiralityConstraint(
  const ChiralityConstraint& chiralityConstraint,
  const Eigen::MatrixXd& positions
) {
  // V = (1 - 4) * [ (2 - 4) x (3 - 4) ]
  return (
    (
      positions.col(std::get<0>(chiralityConstraint))
      - positions.col(std::get<3>(chiralityConstraint))
    ).dot(
      (
        positions.col(std::get<1>(chiralityConstraint))
        - positions.col(std::get<3>(chiralityConstraint))
      ).cross(
        positions.col(std::get<2>(chiralityConstraint))
        - positions.col(std::get<3>(chiralityConstraint))
      )
    )
  );
}

bool moreThanHalfIncorrect(
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

    if( // can this be simplified?
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

template<typename T>
class DGRefinementProblem : public cppoptlib::Problem<T> {
private:
  const std::vector<ChiralityConstraint>& _constraints;

public:
  using typename cppoptlib::Problem<T>::TVector;
  using typename cppoptlib::Problem<T>::THessian;

  DGRefinementProblem(const std::vector<ChiralityConstraint>& constraints) 
    : _constraints(constraints) {}

  T value(const TVector& v); // TODO continue here
  void gradient(const TVector& v, TVector& grad); // TODO continue here

};

bool refine(
  Eigen::MatrixXd& positions,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  DistanceBoundsMatrix& distanceBoundsMatrix
) {
  // embedded is dimensionality x Natoms
  auto numAtoms = positions.cols();
  auto dimensionality = positions.rows();
  assert(dimensionality == 3 || dimensionality == 4);

  /* check if more than half of chirality constraints (which have targets != 0)
   * are incorrect, if so, multiply all z coordinates by -1.
   */
  if(moreThanHalfIncorrect(
    positions,
    chiralityConstraints
  )) {
    positions.row(2) *= -1;
  }
  
}

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip
