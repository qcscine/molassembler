#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "VectorView.h"
#include <cassert>

namespace MoleculeManip {

namespace DistanceGeometry {

DistanceBoundsMatrix::DistanceBoundsMatrix(const unsigned& N) : _N(N) {
  _matrix.resize(_N, _N);
  _matrix.setZero();
  _matrix.triangularView<Eigen::StrictlyUpper>().setConstant(100);

  _initRandomEngine();
}

DistanceBoundsMatrix::DistanceBoundsMatrix(Eigen::MatrixXd matrix) : 
  _matrix(matrix),
  _N(matrix.rows()) {
  assert(matrix.rows() == matrix.cols());

  _initRandomEngine();
}

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

void DistanceBoundsMatrix::triangleInequalitySmooth(
  const SmoothingAlgorithm& algorithmChoice
) {

  if(algorithmChoice == SmoothingAlgorithm::Naive) {
    double ijUpper, ijLower, ikjUpper, ikjLower;

    for(AtomIndexType i = 0; i < _N; i++) {
      for(AtomIndexType j = i + 1; j < _N; j++) {
        ijUpper = upperBound(i, j);
        ijLower = lowerBound(i, j);

        for(AtomIndexType k = j + 1; k < _N; k++) {
          ikjUpper = upperBound(i, k) + upperBound(j, k);
          ikjLower = lowerBound(i, k) + lowerBound(j, k);

          if(ikjUpper < ijUpper) {
            upperBound(i, j) = ikjUpper;
            ijUpper = ikjUpper;
          }

          if(ikjLower > ijLower) {
            lowerBound(i, j) = ikjLower;
            ijLower = ikjLower;
          }
        }
      }
    }
  } else if(algorithmChoice == SmoothingAlgorithm::Custom) {
    /* 1. reformat upper Triangle to i, j, double tuples
     * 2. re-order list to double DESC
     * 3. pick top tuple. use smallest tuples i - k - j to lower i-j bound
     *
     * The idea being that if we improve the highest values first we will have
     * some sort of advantage. This problem I see with the naive implementation
     * is that triangle inequality bounds smoothing since information
     * propagates (at worst) a bond per full run (I think).  Predicated on the
     * idea that information is gained at improvements of large gaps, doing
     * these first gives at least two-bond propagation per run. But this is
     * largely hyperbole, I have no proof. Hence the implementation for
     * testing.
     */

    using BoundTupleType = std::tuple<
      AtomIndexType,
      AtomIndexType,
      double
    >;

    std::vector<BoundTupleType> upperBounds;

    for(AtomIndexType i = 0; i < _N; i++) {
      for(AtomIndexType j = i + 1; j < _N; j++) {
        upperBounds.emplace_back(
          i,
          j,
          upperBound(i, j)
        );
      }
    }

    VectorView<BoundTupleType> sortedView(upperBounds);
    sortedView.sort([](const BoundTupleType& a, const BoundTupleType& b) {
      return std::get<2>(a) > std::get<2>(b);
    });

    VectorView<BoundTupleType> filteredView(upperBounds);

    AtomIndexType i, j;
    double bound;
    for(const BoundTupleType& boundToImprove: sortedView) {
      std::tie(i, j, bound) = boundToImprove;
      
      // filter upperBounds by keeping only tuples containing i OR j
      filteredView.filter([&i, &j](const BoundTupleType& a) -> bool {
        return !(
          (
            std::get<0>(a) == i 
            && std::get<1>(a) != j
          ) || (
            std::get<0>(a) == j
            && std::get<1>(a) != i
          )
        );
      });

      // try all i - k - j  partial sums to see if they're smaller than i - j
      for(AtomIndexType k = 0; k < _N; k++) {
        if(k == i || k == j) continue;

        double sumPartial = 0;
        unsigned found = 0;

        // traverse the tuple list for i-k and j-k
        for(const auto& potentialPartialTuple : filteredView) {
          if( (
              std::get<0>(potentialPartialTuple) == i
              && std::get<1>(potentialPartialTuple) == k
            ) || (
              std::get<0>(potentialPartialTuple) == k
              && std::get<1>(potentialPartialTuple) == i
            ) || (
              std::get<0>(potentialPartialTuple) == j
              && std::get<1>(potentialPartialTuple) == k
            ) || (
              std::get<0>(potentialPartialTuple) == k
              && std::get<1>(potentialPartialTuple) == j
            ) 
          ) {
            sumPartial += std::get<2>(potentialPartialTuple);
            found += 1;
            if(found == 2) break;
          }
        }

        // improves upper bound if smaller!
        if(sumPartial < bound) upperBound(i, j) = sumPartial;
      }

      filteredView.resetFilters(false);
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

  /* Learned some important points from papers:
   * - Going from end to end uniformly negatively affects conformational 
   *   sampling. It is preferable to traverse the list of atoms at random.
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
        lowerBound(i, j),
        upperBound(i, j)
      );

      upperTriangle(
        std::min(i, j),
        std::max(i, j)
      ) = uniformDistribution(_randomEngine);

      /* re-smooth bounds matrix with new information if full metrization is
       * chosen
       */
    }

  }

  return distances;
}

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip
