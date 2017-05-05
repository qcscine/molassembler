#include "DistanceGeometry/DistanceBoundsMatrix.h"
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

DistanceBoundsMatrix::DistanceBoundsMatrix(const Eigen::MatrixXd& matrix) : 
  _boundsMatrix(matrix),
  _N(matrix.rows()) {
  assert(matrix.rows() == matrix.cols());

  _initRandomEngine();
}

/* Private members */
void DistanceBoundsMatrix::_initRandomEngine() {
  std::vector<unsigned> _seeds;

#ifdef NDEBUG
  std::random_device randomDevice;
  for(unsigned n = 0; n < 5; n++) _seeds.emplace_back(randomDevice());
#else 
  _seeds.emplace_back(2721813754);
#endif

  std::seed_seq _seedSequence(_seeds.begin(), _seeds.end());
  _randomEngine.seed(_seedSequence);
}

/* Modifiers */
void DistanceBoundsMatrix::smooth() {
  _boundsMatrix.smooth();
}

bool DistanceBoundsMatrix::setLowerBound(
  const unsigned& i,
  const unsigned& j,
  const double& newLowerBound
) {
  if(
    _boundsMatrix.lowerBound(i, j) < newLowerBound
    && newLowerBound < _boundsMatrix.upperBound(i, j)
  ) {
    _boundsMatrix.lowerBound(i, j) = newLowerBound;
    return true;
  } 

  return false;
}

bool DistanceBoundsMatrix::setUpperBound(
  const unsigned& i,
  const unsigned& j,
  const double& newUpperBound
) {
  if(
    _boundsMatrix.upperBound(i, j) > newUpperBound
    && newUpperBound > _boundsMatrix.lowerBound(i, j)
  ) {
    _boundsMatrix.upperBound(i, j) = newUpperBound;
    return true;
  }

  return false;
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
      ) { // skip on-diagonal and already-chosen elements
        continue; 
      }

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
      if(lowerBound(i, j) > upperBound(i, j)) {
        count += 1;
      }
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

BoundsMatrix DistanceBoundsMatrix::makeSquaredBoundsMatrix() const {
  BoundsMatrix copy = _boundsMatrix;

  for(unsigned i = 0; i < _N; i++) {
    for(unsigned j = i + 1; j < _N; j++) {
      copy.upperBound(i, j) *= copy.upperBound(i, j);
      copy.lowerBound(i, j) *= copy.lowerBound(i, j);
    }
  }

  return copy;
}

double DistanceBoundsMatrix::upperBound(
  const unsigned& i,
  const unsigned& j
) const {
  return _boundsMatrix.upperBound(i, j);
}

} // namespace DistanceGeometry

} // namespace MoleculeManip
