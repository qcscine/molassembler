#ifndef INCLUDE_DG_GENERATE_CONFORMATION_H
#define INCLUDE_DG_GENERATE_CONFORMATION_H

#include "Molecule.h"
#include "Delib/PositionCollection.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/DistanceGeometry.h"

#include <vector>
#include <Eigen/Core>

namespace MoleculeManip {

namespace DistanceGeometry {

namespace detail {

Eigen::Vector3d getPos(
  const Eigen::MatrixXd& positions,
  const AtomIndexType& index
);

double evaluateChiralityConstraint(
  const ChiralityConstraint& chiralityConstraint,
  const Eigen::MatrixXd& positions
);

bool moreThanHalfChiralityConstraintsIncorrect(
  const Eigen::MatrixXd& positions,
  const std::vector<ChiralityConstraint>& chiralityConstraints
);

Delib::PositionCollection convertToPositionCollection(
  const Eigen::VectorXd& vectorizedPositions,
  const EmbeddingOption& embedding
);

template<int segmentSize>
auto getEigen(
  const Eigen::VectorXd& v,
  const unsigned& index, 
  const unsigned& dimensionality
) {
  /* Return a fixed-size const reference to a part of the vector
   *
   * Interesting tidbit to the syntax below:
   * If you write: return v.segment<3>(3 * index);
   *             member function -^^^- integer 3
   *                               |
   *                          operator <
   * 
   * so to disambiguate, you write template before the name of the member
   * function.
   */
  return v.template segment<segmentSize>(dimensionality * index);
}

/* Functor to mimic a HOF that takes a distance Matrix as partial application
 * Using auto makePrototypePropagator(Matrix) {
 *   return [&matrix]()(Prototype) -> ChiralityConstraint {
 *     ...
 *   };
 * }
 * made it impossible to separate header and implementation, as the header is 
 * unusable since auto cannot be deduced.
 */
class PrototypePropagator {
private:
  const Eigen::MatrixXd& distancesMatrix;

public:
  PrototypePropagator(const Eigen::MatrixXd& distancesMatrix);
  ChiralityConstraint operator () (
    const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
  );
};

std::list<Delib::PositionCollection> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
);

} // eo namespace detail

struct MoleculeDGInformation {
  DistanceBoundsMatrix distanceBounds;
  std::vector<
    Stereocenters::Stereocenter::ChiralityConstraintPrototype
  > chiralityConstraintPrototypes;

  MoleculeDGInformation(const unsigned& N);
};

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule
);

// TODO move to namespace detail

// Public functions
std::list<Delib::PositionCollection> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization = MetrizationOption::off,
  const EmbeddingOption& embedding = EmbeddingOption::threeDimensional
);

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization = MetrizationOption::off,
  const EmbeddingOption& embedding = EmbeddingOption::threeDimensional
);

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
