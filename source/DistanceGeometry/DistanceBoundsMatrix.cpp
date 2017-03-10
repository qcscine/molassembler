#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "VectorView.h"
#include <cassert>

namespace MoleculeManip {

namespace DistanceGeometry {

/* Constructors */
DistanceBoundsMatrix::DistanceBoundsMatrix(const unsigned& N) 
  : _boundsMatrix(N),
    _N(N) {

  _boundsMatrix.matrix.triangularView<Eigen::StrictlyUpper>().setConstant(100);
  _initRandomEngine();
}

DistanceBoundsMatrix::DistanceBoundsMatrix(Eigen::MatrixXd matrix) : 
  _boundsMatrix(matrix),
  _N(matrix.rows()) {
  assert(matrix.rows() == matrix.cols());

  _initRandomEngine();
}

/* Private members */
void DistanceBoundsMatrix::_initRandomEngine() {

#ifdef NDEBUG
  std::random_device randomDevice;
  for(unsigned n = 0; n < 5; n++) _seeds.emplace_back(randomDevice());
#else 
  _seeds.emplace_back(2721813754);
#endif

  _seedSequence = std::seed_seq(_seeds.begin(), _seeds.end());
  _randomEngine.seed(_seedSequence);
}

/* Modifiers */
double& DistanceBoundsMatrix::lowerBound(
  const unsigned& i,
  const unsigned& j
) {
  return _boundsMatrix.lowerBound(i, j);
}

void DistanceBoundsMatrix::processDistanceConstraints(
  const std::vector<DistanceConstraint>& constraints
) {
  AtomIndexType i, j;
  double lower, upper;

  for(const auto& constraint : constraints) {
    std::tie(i, j, lower, upper) = constraint;
    /*std::cout << "(" << i << ", " << j << "): [" << lower << ", " << upper << "]"
      << ", currently [" << lowerBound(i, j) << ", " << upperBound(i, j) << "]" 
      << std::endl;*/

    assert(i != j);

    // does applying the constraint reduce slack?
    if(
      upperBound(i, j) > upper // lower constraint
      && upper > lowerBound(i, j) // and it's bigger than the lower bound
    ) {
      upperBound(i, j) = upper;
    }

    if(
      lowerBound(i, j) < lower
      && lower < upperBound(i, j) 
    ) {
      lowerBound(i, j) = lower;
    }
  }
}

double& DistanceBoundsMatrix::upperBound(
  const unsigned& i,
  const unsigned& j
) {
  return _boundsMatrix.upperBound(i, j);
}

void DistanceBoundsMatrix::smooth() {
  _boundsMatrix.smooth();
}

/* Information */
Eigen::MatrixXd DistanceBoundsMatrix::generateDistanceMatrix(
  const MetrizationOption& metrization
) const {
  /* Copy out distanceBounds matrix: Metrization alters bounds to ensure 
   * triangle inequality consistency.
   */
  auto boundsCopy = _boundsMatrix;

  Eigen::MatrixXd distances;
  distances.resize(_N, _N);
  distances.setZero();

  auto upperTriangle = distances.triangularView<Eigen::StrictlyUpper>();

  /* Learned an important point from a Havel paper:
   * Going from end to end negatively affects conformational sampling. It is
   * preferable to traverse the list of atoms at random.
   */

  std::vector<AtomIndexType>  indices(_N);
  std::iota(
    indices.begin(),
    indices.end(),
    0
  );

  std::shuffle(
    indices.begin(),
    indices.end(),
    _randomEngine
  );

  // TODO Using Metrization like this is O(N^5)! 
  for(AtomIndexType idx = 0; idx < _N; idx++) {
    AtomIndexType i = indices.at(idx);
    for(AtomIndexType j = 0; j < _N; j++) {
      if(
        i == j
        || upperTriangle(
          std::min(i, j),
          std::max(i, j)
        ) > 0
      ) continue; // skip on-diagonal and already-chosen elements

      std::uniform_real_distribution<> uniformDistribution(
        boundsCopy.lowerBound(i, j),
        boundsCopy.upperBound(i, j)
      );

      double chosenDistance = uniformDistribution(_randomEngine);

      upperTriangle(
        std::min(i, j),
        std::max(i, j)
      ) = chosenDistance;

      // Full metrization is somewhat naive
      if(metrization == MetrizationOption::full) {
        // Update bounds matrix with chosen value
        boundsCopy.lowerBound(i, j) = chosenDistance;
        boundsCopy.upperBound(i, j) = chosenDistance;

        // Re-smooth the bounds matrix
        boundsCopy.smooth();
      }
    }

  }

  return distances;
}

const Eigen::MatrixXd& DistanceBoundsMatrix::access() const {
  return _boundsMatrix.matrix;
}

unsigned DistanceBoundsMatrix::boundInconsistencies() const {
  unsigned count = 0;

  for(unsigned i = 0; i < _N; i++) {
    for(unsigned j = i + 1; j < _N; j++) {
      if(lowerBound(i, j) > upperBound(i, j)) count++;
    }
  }

  return count;
}

double DistanceBoundsMatrix::lowerBound(
  const unsigned& i,
  const unsigned& j
) const {
  return _boundsMatrix.lowerBound(i, j);
}

double DistanceBoundsMatrix::upperBound(
  const unsigned& i,
  const unsigned& j
) const {
  return _boundsMatrix.upperBound(i, j);
}

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip
