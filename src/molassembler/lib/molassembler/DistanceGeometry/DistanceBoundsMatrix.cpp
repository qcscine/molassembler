#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/Error.h"
#include <cassert>

namespace molassembler {

namespace DistanceGeometry {

constexpr double DistanceBoundsMatrix::defaultLower;
constexpr double DistanceBoundsMatrix::defaultUpper;

DistanceBoundsMatrix::DistanceBoundsMatrix() = default;

DistanceBoundsMatrix::DistanceBoundsMatrix(const unsigned& N) {
  _matrix.resize(N, N);
  _matrix.triangularView<Eigen::Lower>().setConstant(defaultLower);
  _matrix.triangularView<Eigen::StrictlyUpper>().setConstant(defaultUpper);
}

DistanceBoundsMatrix::DistanceBoundsMatrix(Eigen::MatrixXd matrix) : _matrix {std::move(matrix)} {}

bool DistanceBoundsMatrix::setUpperBound(const AtomIndexType i, const AtomIndexType j, const double newUpperBound) {
  if(
    upperBound(i, j) >= newUpperBound
    && newUpperBound > lowerBound(i, j)
  ) {
    _upperBound(i, j) = newUpperBound;
    return true;
  }

  return false;
}

bool DistanceBoundsMatrix::setLowerBound(const AtomIndexType i, const AtomIndexType j, const double newLowerBound) {
  if(
    lowerBound(i, j) <= newLowerBound
    && newLowerBound < upperBound(i, j)
  ) {
    _lowerBound(i, j) = newLowerBound;
    return true;
  }

  return false;
}

void DistanceBoundsMatrix::smooth(Eigen::MatrixXd& matrix) {
  /* Floyd's algorithm: O(N³)
   * Could be refactored slightly that when something is changed, these loops
   * exit and identical loops without the bool setter are run (minimal speed gains)
   */
  const unsigned N = matrix.cols();

  for(AtomIndexType k = 0; k < N; ++k) {
    for(AtomIndexType i = 0; i < N - 1; ++i) {
      for(AtomIndexType j = i + 1; j < N; ++j) {
        if(upperBound(matrix, i, j) > upperBound(matrix, i, k) + upperBound(matrix, k, j)) {
          upperBound(matrix, i, j) = upperBound(matrix, i, k) + upperBound(matrix, k, j);
        }

        if(lowerBound(matrix, i, j) < lowerBound(matrix, i, k) - upperBound(matrix, k, j)) {
          lowerBound(matrix, i, j) = lowerBound(matrix, i, k) - upperBound(matrix, k, j);
        } else if(lowerBound(matrix, i, j) < lowerBound(matrix, j, k) - upperBound(matrix, k, i)) {
          lowerBound(matrix, i, j) = lowerBound(matrix, j, k) - upperBound(matrix, k, i);
        }

        assert(lowerBound(matrix, i, j) <= upperBound(matrix, i, j));
      }
    }
  }
}

void DistanceBoundsMatrix::smooth() {
  smooth(_matrix);
}

unsigned DistanceBoundsMatrix::boundInconsistencies() const {
  unsigned count = 0;
  const unsigned N = _matrix.cols();

  for(unsigned i = 0; i < N - 1; i++) {
    for(unsigned j = i + 1; j < N; j++) {
      if(lowerBound(i, j) > upperBound(i, j)) {
        count += 1;
      }
    }
  }

  return count;
}

const Eigen::MatrixXd& DistanceBoundsMatrix::access() const {
  return _matrix;
}

outcome::result<Eigen::MatrixXd> DistanceBoundsMatrix::makeDistanceMatrix(Partiality partiality) const noexcept {
  auto matrixCopy = _matrix;

  const unsigned N = _matrix.cols();

  std::vector<AtomIndexType> indices(N);
  std::iota(
    indices.begin(),
    indices.end(),
    0
  );

  temple::random.shuffle(indices);

  std::vector<AtomIndexType>::const_iterator separator;

  if(partiality == Partiality::FourAtom) {
    separator = indices.cbegin() + std::min(N, 4u);
  } else if(partiality == Partiality::TenPercent) {
    separator = indices.cbegin() + std::min(N, static_cast<unsigned>(0.1 * N));
  } else { // All
    separator = indices.cend();
  }

  for(auto iter = indices.cbegin(); iter != separator; ++iter) {
    const AtomIndexType i = *iter;
    for(AtomIndexType j = 0; j < N; j++) {
      if(
        i == j
        || matrixCopy(i, j) == matrixCopy(j, i)
      ) { // skip on-diagonal and already-chosen elements
        continue;
      }

      if(upperBound(matrixCopy, i, j) < lowerBound(matrixCopy, i, j)) {
        return DGError::GraphImpossible;
      }

      double chosenDistance = temple::random.getSingle<double>(
        lowerBound(matrixCopy, i, j),
        upperBound(matrixCopy, i, j)
      );

      lowerBound(matrixCopy, i, j) = chosenDistance;
      upperBound(matrixCopy, i, j) = chosenDistance;

      smooth(matrixCopy);
    }
  }

  for(auto iter = separator; iter != indices.cend(); ++iter) {
    const AtomIndexType i = *iter;
    for(AtomIndexType j = 0; j < N; j++) {
      if(
        i == j
        || matrixCopy(i, j) == matrixCopy(j, i)
      ) { // skip on-diagonal and already-chosen elements
        continue;
      }

      /* Interval reversal is no longer a failure criterion, we have already
       * thrown caution to the wind by no longer smoothing the matrix.
       * Nevertheless, to avoid UB, it is still necessary to properly order
       * the parameters
       */
      double chosenDistance = temple::random.getSingle<double>(
        std::min(
          lowerBound(matrixCopy, i, j),
          upperBound(matrixCopy, i, j)
        ),
        std::max(
          lowerBound(matrixCopy, i, j),
          upperBound(matrixCopy, i, j)
        )
      );

      lowerBound(matrixCopy, i, j) = chosenDistance;
      upperBound(matrixCopy, i, j) = chosenDistance;
    }
  }

  return matrixCopy;
}

Eigen::MatrixXd DistanceBoundsMatrix::makeSquaredBoundsMatrix() const {
  return _matrix.cwiseProduct(_matrix);
}

unsigned DistanceBoundsMatrix::N() const {
  return _matrix.cols();
}

} // namespace DistanceGeometry

} // namespace molassembler
