/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"

#include "molassembler/DistanceGeometry/Error.h"
#include "molassembler/Options.h"

#include "temple/Random.h"

#include <cassert>
#include <algorithm>

namespace Scine {

namespace molassembler {

namespace DistanceGeometry {

constexpr double DistanceBoundsMatrix::defaultLower;
constexpr double DistanceBoundsMatrix::defaultUpper;

DistanceBoundsMatrix::DistanceBoundsMatrix() = default;

DistanceBoundsMatrix::DistanceBoundsMatrix(const long unsigned N) {
  _matrix.resize(N, N);
  _matrix.triangularView<Eigen::Lower>().setConstant(defaultLower);
  _matrix.triangularView<Eigen::StrictlyUpper>().setConstant(defaultUpper);
}

DistanceBoundsMatrix::DistanceBoundsMatrix(Eigen::MatrixXd matrix) : _matrix {std::move(matrix)} {}

bool DistanceBoundsMatrix::setUpperBound(const AtomIndex i, const AtomIndex j, const double newUpperBound) {
  if(
    upperBound(i, j) >= newUpperBound
    && newUpperBound > lowerBound(i, j)
  ) {
    _upperBound(i, j) = newUpperBound;
    return true;
  }

  return false;
}

bool DistanceBoundsMatrix::setLowerBound(const AtomIndex i, const AtomIndex j, const double newLowerBound) {
  if(
    lowerBound(i, j) <= newLowerBound
    && newLowerBound < upperBound(i, j)
  ) {
    _lowerBound(i, j) = newLowerBound;
    return true;
  }

  return false;
}

void DistanceBoundsMatrix::smooth(Eigen::Ref<Eigen::MatrixXd> matrix) {
  /* Floyd's algorithm: O(N³) */
  const unsigned N = matrix.cols();

  for(AtomIndex k = 0; k < N; ++k) {
    for(AtomIndex i = 0; i < N - 1; ++i) {
      /* Single-branch references to lower and upper ik parts */
      auto upperLowerIK = (i < k
        ? std::pair<double&, double&>(matrix(i, k), matrix(k, i))
        : std::pair<double&, double&>(matrix(k, i), matrix(i, k))
      );
      double& upperIK = upperLowerIK.first;
      double& lowerIK = upperLowerIK.second;

      if(lowerIK > upperIK) {
        throw std::runtime_error("Triangle smoothing encountered bound inversion");
      }

      for(AtomIndex j = i + 1; j < N; ++j) {
        /* i < j is known, so lower(matrix, i, j) is always matrix(j, i),
         * and upper(matrix, i, j) always matrix(i, j)
         */
        double& upperIJ = matrix(i, j);
        double& lowerIJ = matrix(j, i);

        /* Single-branch references to lower and upper jk parts */
        auto upperLowerJK = (j < k
          ? std::pair<double&, double&>(matrix(j, k), matrix(k, j))
          : std::pair<double&, double&>(matrix(k, j), matrix(j, k))
        );
        double& upperJK = upperLowerJK.first;
        double& lowerJK = upperLowerJK.second;

        /* Actual algorithm */
        if(upperIJ > upperIK + upperJK) {
          upperIJ = upperIK + upperJK;
        }

        if(lowerIJ < lowerIK - upperJK) {
          lowerIJ = lowerIK - upperJK;
        } else if(lowerIJ < lowerJK - upperIK) {
          lowerIJ = lowerJK - upperIK;
        }

        // Safety
        if(lowerIJ > upperIJ) {
          throw std::runtime_error("Triangle smoothing encountered bound inversion");
        }
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

outcome::result<Eigen::MatrixXd> DistanceBoundsMatrix::makeDistanceMatrix(random::Engine& engine, Partiality partiality) const noexcept {
  auto matrixCopy = _matrix;

  const unsigned N = _matrix.cols();

  std::vector<AtomIndex> indices(N);
  std::iota(
    indices.begin(),
    indices.end(),
    0
  );

  temple::random::shuffle(indices, engine);

  std::vector<AtomIndex>::const_iterator separator;

  if(partiality == Partiality::FourAtom) {
    separator = indices.cbegin() + std::min(N, 4u);
  } else if(partiality == Partiality::TenPercent) {
    separator = indices.cbegin() + std::min(N, static_cast<unsigned>(0.1 * N));
  } else { // All
    separator = indices.cend();
  }

  // Up to the separator decided by partiality
  for(auto iter = indices.cbegin(); iter != separator; ++iter) {
    const AtomIndex i = *iter;
    for(AtomIndex j = 0; j < N; j++) {
      if(
        i == j
        || matrixCopy(i, j) == matrixCopy(j, i)
      ) { // skip on-diagonal and already-chosen elements
        continue;
      }

      if(upperBound(matrixCopy, i, j) < lowerBound(matrixCopy, i, j)) {
        return DgError::GraphImpossible;
      }

      double chosenDistance = temple::random::getSingle<double>(
        lowerBound(matrixCopy, i, j),
        upperBound(matrixCopy, i, j),
        engine
      );

      lowerBound(matrixCopy, i, j) = chosenDistance;
      upperBound(matrixCopy, i, j) = chosenDistance;

      smooth(matrixCopy);
    }
  }

  for(auto iter = separator; iter != indices.cend(); ++iter) {
    const AtomIndex i = *iter;
    for(AtomIndex j = 0; j < N; j++) {
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
       *
       * TODO this is an atrocious access pattern (see impl of lowerBound)
       */
      double chosenDistance = temple::random::getSingle<double>(
        std::min(
          lowerBound(matrixCopy, i, j),
          upperBound(matrixCopy, i, j)
        ),
        std::max(
          lowerBound(matrixCopy, i, j),
          upperBound(matrixCopy, i, j)
        ),
        engine
      );

      lowerBound(matrixCopy, i, j) = chosenDistance;
      upperBound(matrixCopy, i, j) = chosenDistance;
    }
  }

  return matrixCopy;
}

Eigen::MatrixXd DistanceBoundsMatrix::makeSquaredBoundsMatrix() const {
  return _matrix.array().square();
}

unsigned DistanceBoundsMatrix::N() const {
  return _matrix.cols();
}

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine
