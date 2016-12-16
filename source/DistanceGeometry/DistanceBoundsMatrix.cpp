#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include <cassert>

namespace MoleculeManip {

namespace DistanceGeometry {

DistanceBoundsMatrix::DistanceBoundsMatrix(const unsigned& N) : _N(N) {
  _matrix.resize(_N, _N);
  _matrix.setZero();
  _matrix.triangularView<Eigen::StrictlyUpper>().setConstant(100);

#ifdef NDEBUG
  std::random_device randomDevice;
  for(unsigned n = 0; n < 5; n++) _seeds.emplace_back(randomDevice());
#else 
  _seeds.emplace_back(2721813754);
#endif
  _seedSequence = std::seed_seq(_seeds.begin(), _seeds.end());
  _randomEngine.seed(_seedSequence);
}

decltype(DistanceBoundsMatrix::_matrix(1, 2))& DistanceBoundsMatrix::upperBound(
  const unsigned& i,
  const unsigned& j
) {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}

decltype(DistanceBoundsMatrix::_matrix(1, 2))& DistanceBoundsMatrix::lowerBound(
  const unsigned& i,
  const unsigned& j
) {
  return _matrix(
    std::max(i, j),
    std::min(i, j)
  );
}

double DistanceBoundsMatrix::upperBound(
  const unsigned& i,
  const unsigned& j
) const {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}

double DistanceBoundsMatrix::lowerBound(
  const unsigned& i,
  const unsigned& j
) const {
  return _matrix(
    std::max(i, j),
    std::min(i, j)
  );
}

void DistanceBoundsMatrix::processDistanceConstraints(
  const std::vector<DistanceConstraint>& constraints
) {
  for(const auto& constraint : constraints) {
    AtomIndexType i, j;
    double lower, upper;

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

// TODO alter behavior due to metrizationOption!
Eigen::MatrixXd DistanceBoundsMatrix::generateDistanceMatrix(
  const MetrizationOption& metrization
) const { // cannot be const since it alters the randomEngine state!

  Eigen::MatrixXd distances;
  distances.resize(_N, _N);
  distances.setZero();

  auto upperTriangle = distances.triangularView<Eigen::StrictlyUpper>();

  for(unsigned i = 0; i < _N; i++) {
    for(unsigned j = i + 1; j < _N; j++) {
      std::uniform_real_distribution<> uniformDistribution(
        lowerBound(i, j),
        upperBound(i, j)
      );
      upperTriangle(i, j) = uniformDistribution(_randomEngine);
    }
  }

  return distances;
}

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip
