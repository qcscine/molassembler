/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "chemical_symmetries/Recognition.h"

#include "chemical_symmetries/CoordinateSystemTransformation.h"
#include "chemical_symmetries/Diophantine.h"
#include "chemical_symmetries/Partitioner.h"

#include "boost/optional.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include "temple/Functional.h"
#include "temple/Optimization/Common.h"
#include "temple/Optimization/NelderMead.h"
#include "temple/Optimization/TrustRegion.h"

#include "temple/Stringify.h"
#include <iostream>
#include <fstream>
#include <random>

/* TODO
 * - standardizeTop needs to orient asymmetric tops along z with the highest
 *   order cn axis found
 * - Perf change if i-j loop break in greedy variants to permutational
 *   calculations is a continue instead (both break and continue are viable as
 *   long as the swap-back is skipped)
 * - Greedy variants are close but not quite there. Perhaps 2x or start from
 *   random permutation works better? Maybe it's also not worth it in the
 *   diophantine functions, only in the all elements permutation. Time to
 *   benchmark, I think.
 * - Allowing npGroupings for l = G / np = 1 has tanked performance heavily
 *   Maybe need to be smarter...
 */

namespace Scine {
namespace Symmetry {

namespace detail {

Eigen::MatrixXd dropColumn(const Eigen::MatrixXd& matrix, unsigned colToRemove) {
  const unsigned c = matrix.cols();
  const unsigned leftColCount = colToRemove - 1;
  const unsigned rightColCount = c - colToRemove;

  Eigen::MatrixXd copy(matrix.rows(), c - 1);
  copy.leftCols(leftColCount) = matrix.leftCols(leftColCount);
  copy.rightCols(rightColCount) = matrix.rightCols(rightColCount);

  return copy;
}

/*! @brief Normalize positions for continuous symmetry measure analysis
 *
 * Reframes to center of mass frame (although no masses exist) and rescales
 * vectors so that the maximum distance is 1)
 */
PositionCollection normalize(const PositionCollection& positions) {
  const unsigned N = positions.cols();

  // Translate the origin to the average position
  const Eigen::Vector3d center = positions.rowwise().sum() / positions.cols();
  PositionCollection transformed = positions.colwise() - center;

  // Rescale all distances so that the longest is a unit vector
  const double longestDistance = std::sqrt(positions.colwise().squaredNorm().maxCoeff());
  for(unsigned i = 0; i < N; ++i) {
    transformed.col(i) /= longestDistance;
  }

  // Drop any vectors that are very close to the center of mass
  for(int i = 0; i < transformed.cols(); ++i) {
    if(transformed.col(i).norm() < 1e-3) {
      const int r = transformed.rows();
      const int c = transformed.cols() - 1;
      if(i < c) {
        transformed.block(0, i, r, c - i) = transformed.rightCols(c - i);
      }
      transformed.conservativeResize(r, c);
      --i;
    }
  }

  assert(transformed.cols() >= 2);

  return transformed;
}

//! Determine degeneracy of intertial moments
unsigned degeneracy(const Eigen::Vector3d& inertialMoments) {
  constexpr double degeneracyEpsilon = 0.05;
  unsigned mdeg = 0;
  if(
    std::fabs(
      (inertialMoments(2) - inertialMoments(1)) / inertialMoments(2)
    ) <= degeneracyEpsilon
  ) {
    mdeg += 1;
  }
  if(
    std::fabs(
      (inertialMoments(1) - inertialMoments(0)) / inertialMoments(2)
    ) <= degeneracyEpsilon
  ) {
    mdeg += 2;
  }

  return 1 + (mdeg + 1) / 2;
}

//! Determine point group of an approximately linear system
PointGroup linear(const PositionCollection& normalizedPositions) {
  const unsigned P = normalizedPositions.cols();
  for(unsigned i = 0; i < P / 2; ++i) {
    const unsigned opposite = P - i - 1;
    if(!normalizedPositions.col(i).isApprox(-normalizedPositions.col(opposite), 1e-4)) {
      /* Positions do not have inversion symmetry */
      return PointGroup::Cinfv;
    }
  }

  return PointGroup::Dinfh;
}

} // namespace detail

namespace csm {

struct IdentityMap {
  PURITY_STRONG inline unsigned at(const unsigned i) const noexcept {
    return i;
  }
};

double calculateCSM(
  const PositionCollection& normalizedPositions,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  const std::vector<unsigned>& particles
) {
  const unsigned p = particles.size();

  double value = 0;

  /* Fold points and average */
  Eigen::Vector3d averagePoint = Eigen::Vector3d::Zero();
  for(unsigned i = 0; i < p; ++i) {
    averagePoint += foldMatrices.block<3, 3>(0, 3 * i) * normalizedPositions.col(particles.at(i));
  }
  averagePoint /= p;

  /* Unfold points */
  for(unsigned i = 0; i < p; ++i) {
    value += (
      unfoldMatrices.block<3, 3>(0, 3 * i) * averagePoint
      - normalizedPositions.col(particles.at(i))
    ).squaredNorm();
  }

  value *= 100.0 / p;
  return value;
}

double calculateCSM(
  const PositionCollection& normalizedPositions,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  const std::vector<unsigned>& particles,
  const elements::ElementGrouping& elementGrouping
) {
  const unsigned p = particles.size();
  const unsigned l = elementGrouping.groups.front().size();
  assert(p * l == unfoldMatrices.cols() / 3);

  double value = 0;

  /* Fold points and average */
  Eigen::Vector3d averagePoint = Eigen::Vector3d::Zero();
  for(unsigned i = 0; i < p; ++i) {
    const auto& elements = elementGrouping.groups.at(i);
    for(unsigned j = 0; j < l; ++j) {
      averagePoint += foldMatrices.block<3, 3>(0, 3 * elements.at(j)) * normalizedPositions.col(particles.at(i));
    }
  }
  averagePoint /= (p * l);

  /* Unfold points */
  for(unsigned i = 0; i < p; ++i) {
    const auto& elements = elementGrouping.groups.at(i);
    for(unsigned j = 0; j < l; ++j) {
      value += (
        unfoldMatrices.block<3, 3>(0, 3 * elements.at(j)) * averagePoint
        - normalizedPositions.col(particles.at(i))
      ).squaredNorm();
    }
  }

  value *= 100.0 / (p * l);
  return value;
}


double allSymmetryElements(
  const PositionCollection& normalizedPositions,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  std::vector<unsigned> particleIndices
) {
  const unsigned P = particleIndices.size();
  const unsigned G = unfoldMatrices.cols() / 3;
  assert(P == G);

  double value = 1000;

  do {
    double permutationCSM = calculateCSM(
      normalizedPositions,
      unfoldMatrices,
      foldMatrices,
      particleIndices
    );

    value = std::min(value, permutationCSM);
  } while(
    std::next_permutation(
      std::begin(particleIndices),
      std::end(particleIndices)
    )
  );

  return value;
}

double greedyAllSymmetryElements(
  const PositionCollection& normalizedPositions,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  std::vector<unsigned> particleIndices
) {
  const unsigned P = particleIndices.size();
  const unsigned G = unfoldMatrices.cols() / 3;
  assert(P == G);

  /* Each search for a better adjacent permutation is O(N^2) worst case
   * The number of searches for a better adjacent permutation is at least on
   * the order of O(N), so let's roughly estimate this as O(N^3)
   */
  bool foundBetterAdjacentPermutation = false;
  double currentCSM = calculateCSM(
    normalizedPositions,
    unfoldMatrices,
    foldMatrices,
    particleIndices
  );
  do {
    foundBetterAdjacentPermutation = false;
    for(unsigned i = 0; i < G && !foundBetterAdjacentPermutation; ++i) {
      for(unsigned j = i + 1; j < G; ++j) {
        std::swap(particleIndices.at(i), particleIndices.at(j));

        double adjacentCSM = calculateCSM(
          normalizedPositions,
          unfoldMatrices,
          foldMatrices,
          particleIndices
        );
        if(adjacentCSM < currentCSM) {
          currentCSM = adjacentCSM;
          foundBetterAdjacentPermutation = true;
          continue;
        }

        std::swap(particleIndices.at(i), particleIndices.at(j));
      }
    }
  } while(foundBetterAdjacentPermutation);

  return currentCSM;
}

double groupedSymmetryElements(
  const PositionCollection& normalizedPositions,
  std::vector<unsigned> particleIndices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  const elements::ElementGrouping& elementGrouping
) {
  /* The number of groups in element grouping must match the number of particles
   * permutated here
   */
  assert(elementGrouping.groups.size() == particleIndices.size());
  assert(std::is_sorted(std::begin(particleIndices), std::end(particleIndices)));

  double value = 1000;

  do {
    double permutationCSM = calculateCSM(
      normalizedPositions,
      unfoldMatrices,
      foldMatrices,
      particleIndices,
      elementGrouping
    );
    value = std::min(value, permutationCSM);
  } while(
    std::next_permutation(
      std::begin(particleIndices),
      std::end(particleIndices)
    )
  );

  return value;
}

double greedyGroupedSymmetryElements(
  const PositionCollection& normalizedPositions,
  std::vector<unsigned> particleIndices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  const elements::ElementGrouping& elementGrouping
) {
  /* The number of groups in element grouping must match the number of particles
   * permutated here
   */
  assert(elementGrouping.groups.size() == particleIndices.size());
  assert(std::is_sorted(std::begin(particleIndices), std::end(particleIndices)));

  const unsigned groupSize = particleIndices.size();

  bool foundBetterAdjacentPermutation = false;
  double currentCSM = calculateCSM(
    normalizedPositions,
    unfoldMatrices,
    foldMatrices,
    particleIndices,
    elementGrouping
  );
  do {
    foundBetterAdjacentPermutation = false;
    for(unsigned i = 0; i < groupSize && !foundBetterAdjacentPermutation; ++i) {
      for(unsigned j = i + 1; j < groupSize; ++j) {
        std::swap(particleIndices.at(i), particleIndices.at(j));

        double adjacentCSM = calculateCSM(
          normalizedPositions,
          unfoldMatrices,
          foldMatrices,
          particleIndices,
          elementGrouping
        );
        if(adjacentCSM < currentCSM) {
          currentCSM = adjacentCSM;
          foundBetterAdjacentPermutation = true;
          continue;
        }

        std::swap(particleIndices.at(i), particleIndices.at(j));
      }
    }
  } while(foundBetterAdjacentPermutation);

  return currentCSM;
}

struct OrientationCSMFunctor {
  using MatrixType = Eigen::Matrix<double, 3, Eigen::Dynamic>;

  const PositionCollection& coordinates;
  const MatrixType unfoldMatrices;
  const MatrixType foldMatrices;
  const std::unordered_map<unsigned, elements::ElementGrouping> npGroups;

  OrientationCSMFunctor(
    const PositionCollection& normalizedPositions,
    const PointGroup group
  ) : coordinates(normalizedPositions),
      unfoldMatrices(makeUnfoldMatrices(elements::symmetryElements(group))),
      foldMatrices(makeFoldMatrices(unfoldMatrices)),
      npGroups(elements::npGroupings(elements::symmetryElements(group)))
  {}

  static MatrixType makeUnfoldMatrices(const elements::ElementsList& elements) {
    MatrixType unfold(3, elements.size() * 3);
    const unsigned G = elements.size();
    for(unsigned i = 0; i < G; ++i) {
      unfold.block<3, 3>(0, 3 * i) = elements.at(i)->matrix();
    }
    return unfold;
  }

  static MatrixType makeFoldMatrices(const MatrixType& unfoldMatrices) {
    MatrixType fold(3, unfoldMatrices.cols());
    const unsigned G = unfoldMatrices.cols() / 3;
    for(unsigned i = 0; i < G; ++i) {
      fold.block<3, 3>(0, 3 * i) = unfoldMatrices.block<3, 3>(0, 3 * i).inverse();
    }
    return fold;
  }

  static Eigen::Matrix3d eulerRotation(const Eigen::VectorXd& parameters) {
    // Three euler angles in radians
    assert(parameters.size() == 3);
    const double alpha = parameters(0);
    const double beta = parameters(1);
    const double gamma = parameters(2);

    // Figure out N and z'
    const Eigen::Vector3d N = Eigen::AngleAxisd(
      alpha,
      Eigen::Vector3d::UnitZ()
    ) * Eigen::Vector3d::UnitX();

    const Eigen::Vector3d zPrime = Eigen::AngleAxisd(
      beta,
      N
    ) * Eigen::Vector3d::UnitZ();

    // Compose the euler rotation matrix
    return (
      Eigen::AngleAxisd(gamma, zPrime.normalized())
      * Eigen::AngleAxisd(beta, N.normalized())
      * Eigen::AngleAxisd(alpha, Eigen::Vector3d::UnitZ())
    ).toRotationMatrix();
  }

  double direct_csm(
    const PositionCollection& positions,
    const std::vector<unsigned>& particleIndices
  ) const {
    assert(particleIndices.size() == static_cast<decltype(particleIndices.size())>(foldMatrices.cols()) / 3);
    return allSymmetryElements(
      positions,
      unfoldMatrices,
      foldMatrices,
      particleIndices
    );
  }

  double diophantine_csm(
    const PositionCollection& positions,
    const std::vector<unsigned>& subdivisionGroupSizes,
    const std::vector<unsigned>& particleIndices
  ) const {
    const unsigned P = particleIndices.size();
    std::vector<unsigned> subdivisionMultipliers;
    if(!diophantine::first_solution(subdivisionMultipliers, subdivisionGroupSizes, P)) {
      throw std::logic_error("Diophantine failure! Couldn't find first solution");
    }

    double value = 1000;
    do {
      /* We have a composition of groups to subdivide our points:
       *
       * E.g.: P = 12
       * subdivisionGroupSizes {4, 3, 2}
       * one solution: subdivisionMultipliers {1, 2, 1}
       *
       * Next we make every possible partition of our points into those
       * respective group sizes. We populate a flat map of point index to
       * group size index according to the group sizes and multipliers:
       *
       * p0 = 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2
       *
       * Of this, there are 12! / (4! 6! 2!) = 13860 permutations.
       */
      std::vector<unsigned> flatGroupMap;
      flatGroupMap.reserve(P);
      for(unsigned i = 0; i < subdivisionMultipliers.size(); ++i) {
        const unsigned groupMultiplier = subdivisionMultipliers.at(i);

        if(groupMultiplier == 0) {
          continue;
        }

        const unsigned groupSize = subdivisionGroupSizes.at(i);
        flatGroupMap.resize(flatGroupMap.size() + groupSize * groupMultiplier, i);
      }
      assert(flatGroupMap.size() == P);

      do {
        /* For each of these permutations, we sub-partition each group using
         * Partitioner. This is necessary to treat multipliers > 1 correctly.
         *
         * In the previous case, where we have a two groups of size three, we
         * have 10 sub-partitions to consider. All other groups with
         * multiplier == 1 have only a single sub-partition.
         */

        /* Collect the indices mapped to groups of equal size */
        const unsigned numSizeGroups = subdivisionGroupSizes.size();
        std::vector<
          std::vector<unsigned>
        > sameSizeIndexGroups(numSizeGroups);
        for(unsigned i = 0; i < P; ++i) {
          sameSizeIndexGroups.at(flatGroupMap.at(i)).push_back(i);
        }

        double partitionCSM = 0;
        for(unsigned i = 0; i < numSizeGroups; ++i) {
          const unsigned groupSize = subdivisionGroupSizes.at(i);
          const unsigned multiplier = subdivisionMultipliers.at(i);
          const auto& sameSizeParticleIndices = sameSizeIndexGroups.at(i);
          const auto& npGroup = npGroups.at(groupSize);

          // Skip groups with multiplier zero
          if(multiplier == 0) {
            continue;
          }

          /* Subpartition csm is minimized over the sub-partitions of
           * groups of equal size
           */
          Partitioner partitioner {multiplier, groupSize};
          double bestSubpartitionCSM = 1000;
          do {
            double subpartitionCSM = 0;
            for(auto&& partitionIndices : partitioner.partitions()) {
              /* Map all the way back to actual particle indices:
               * partition -> same size particles -> particle
               */
              auto partitionParticles = temple::map(
                partitionIndices,
                [&](const unsigned indexOfParticle) -> unsigned {
                  return particleIndices.at(
                    sameSizeParticleIndices.at(indexOfParticle)
                  );
                }
              );

              subpartitionCSM += groupedSymmetryElements(
                positions,
                partitionParticles,
                unfoldMatrices,
                foldMatrices,
                npGroup
              );
            }
            bestSubpartitionCSM = std::min(bestSubpartitionCSM, subpartitionCSM);
          } while(partitioner.next_partition());

          partitionCSM += bestSubpartitionCSM;
        }

        value = std::min(value, partitionCSM);
      } while(std::next_permutation(std::begin(flatGroupMap), std::end(flatGroupMap)));
    } while(diophantine::next_solution(subdivisionMultipliers, subdivisionGroupSizes, P));

    return value;
  }

  double csm(const PositionCollection& positions) const {
    const unsigned P = positions.cols();
    const unsigned G = foldMatrices.cols() / 3;

    if(P == G) {
      return direct_csm(
        positions,
        temple::iota<unsigned>(P)
      );
    }

    if(P == G + 1) {
      double minimalCSM = 1000;

      for(unsigned i = 0; i < P; ++i) {
        std::vector<unsigned> particleIndices;
        particleIndices.reserve(G);
        for(unsigned j = 0; j < i; ++j) {
          particleIndices.push_back(j);
        }
        for(unsigned j = i + 1; j < P; ++j) {
          particleIndices.push_back(j);
        }

        const double oneSymmetrizedToOriginCSM = (
          direct_csm(positions, particleIndices)
          + positions.col(i).squaredNorm() * 100
        );

        minimalCSM = std::min(minimalCSM, oneSymmetrizedToOriginCSM);
      }

      return minimalCSM;
    }

    /* Set up list of subdivisionGroupSizes */
    std::vector<unsigned> subdivisionGroupSizes {};
    /* The sizes of groups that we can subdivide contains the full set only
     * if there are more particles than symmetry elements:
     */
    if(P > G) {
      subdivisionGroupSizes.push_back(G);
    }
    for(const auto& groupMapPair : npGroups) {
      subdivisionGroupSizes.push_back(groupMapPair.first);
    }
    std::sort(
      std::begin(subdivisionGroupSizes),
      std::end(subdivisionGroupSizes),
      std::greater<>()
    );
    if(diophantine::has_solution(subdivisionGroupSizes, P)) {
      return diophantine_csm(
        positions,
        subdivisionGroupSizes,
        temple::iota<unsigned>(P)
      );
    }

    if(P > 2 && diophantine::has_solution(subdivisionGroupSizes, P - 1)) {
      double minimalCSM = 1000;

      for(unsigned i = 0; i < P; ++i) {
        std::vector<unsigned> particleIndices;
        particleIndices.reserve(G);
        for(unsigned j = 0; j < i; ++j) {
          particleIndices.push_back(j);
        }
        for(unsigned j = i + 1; j < P; ++j) {
          particleIndices.push_back(j);
        }

        const double oneSymmetrizedToOriginCSM = (
          diophantine_csm(positions, subdivisionGroupSizes, particleIndices)
          + positions.col(i).squaredNorm() * 100
        );

        minimalCSM = std::min(minimalCSM, oneSymmetrizedToOriginCSM);
      }

      return minimalCSM;
    }

    throw std::logic_error("You shouldn't even instantiate this type if you know that you cannot calculate a CSM");
  }

  double operator() (const Eigen::VectorXd& parameters) const {
    const PositionCollection rotatedCoordinates = eulerRotation(parameters) * coordinates;
    return csm(rotatedCoordinates);
  }
};

boost::optional<double> pointGroup(
  const PositionCollection& normalizedPositions,
  const PointGroup group
) {
  const auto elements = elements::symmetryElements(group);
  const auto npGroups = npGroupings(elements);

  const unsigned G = elements.size();
  const unsigned P = normalizedPositions.cols();

  /* There are conditions on when we can calculate a CSM for this number of
   * particles P and the particular point group with G symmetry elements.
   */

  /* Set up a list of group sizes we can subdivide particles into */
  std::vector<unsigned> subdivisionGroupSizes {};
  /* The sizes of groups that we can subdivide contains the full set only
   * if there are more particles than symmetry elements:
   */
  if(P > G) {
    subdivisionGroupSizes.push_back(G);
  }
  for(const auto& groupMapPair : npGroups) {
    subdivisionGroupSizes.push_back(groupMapPair.first);
  }
  std::sort(
    std::begin(subdivisionGroupSizes),
    std::end(subdivisionGroupSizes),
    std::greater<>()
  );

  if(
    !(P == G)
    && !(P == G + 1)
    && !diophantine::has_solution(subdivisionGroupSizes, P)
    && !(
      P > 2
      && diophantine::has_solution(subdivisionGroupSizes, P - 1)
    )
  ) {
    /* If none of these possibilities are true, then we cannot calculate a CSM
     * for this number of particles and this particular point group.
     */
    return boost::none;
  }

  OrientationCSMFunctor functor {normalizedPositions, group};

  // Set up the initial simplex
  Eigen::MatrixXd simplex(3, 4);
  simplex.col(0) = Eigen::Vector3d::Zero();
  simplex.block<3, 3>(0, 1) = Eigen::Matrix3d::Identity();

  struct NelderMeadChecker {
    bool shouldContinue(unsigned iteration, double lowestValue, double stddev) const {
      return iteration < 1000 && lowestValue > 1e-3 && stddev > 1e-2;
    }
  };

  auto minimizationResult = temple::NelderMead<>::minimize(
    simplex,
    functor,
    NelderMeadChecker {}
  );

  std::cout << "Minimized to " << minimizationResult.value << " in " << minimizationResult.iterations << " iterations. Final simplex:\n " << simplex << "\n";

  return minimizationResult.value;
}

namespace cn {

PositionCollection symmetrize(
  const PositionCollection& normalizedPositions,
  const Eigen::Vector3d& axis,
  const unsigned n,
  const std::vector<
    std::vector<unsigned>
  > groups
) {
  PositionCollection symmetrizedPositions = normalizedPositions;

  for(const auto& group : groups) {
    const double angleRadians = 2 * M_PI / n;
    /* Fold */
    // The first position is unchanged, hence we start with i = 1
    for(unsigned i = 1; i < n; ++i) {
      symmetrizedPositions.col(group.front()) += Eigen::AngleAxisd(i * angleRadians, axis) * normalizedPositions.col(group.at(i));
    }

    // Average the folded positions
    symmetrizedPositions.col(group.front()) /= n;

    /* Unfold */
    for(unsigned i = 1; i < n; ++i) {
      symmetrizedPositions.col(group.at(i)) = Eigen::AngleAxisd(-i * angleRadians, axis) * symmetrizedPositions.col(group.front());
    }
  }

  return symmetrizedPositions;
}

double evaluatePermutation(
  const PositionCollection& normalizedPositions,
  const Eigen::Vector3d& axis,
  const unsigned n,
  const std::vector<unsigned>& groupPermutation
) {
  assert(n >= 2);

  const double angleRadians = 2 * M_PI / n;
  Eigen::Vector3d averagePoint = normalizedPositions.col(groupPermutation.front());

  /* Fold and average */
  averagePoint = normalizedPositions.col(groupPermutation.front());
  for(unsigned i = 1; i < n; ++i) {
    averagePoint.noalias() += Eigen::AngleAxisd(i * angleRadians, axis) * normalizedPositions.col(groupPermutation.at(i));
  }
  averagePoint /= n;

  /* Calculate CSM while unfolding */
  double csm = (
    normalizedPositions.col(groupPermutation.front())
    - averagePoint
  ).squaredNorm();

  for(unsigned i = 1; i < n; ++i) {
    csm += (
      normalizedPositions.col(groupPermutation.at(i))
      - Eigen::AngleAxisd(i * angleRadians, -axis) * averagePoint
    ).squaredNorm();
  }

  return csm;
}

/**
 * @brief Functor for minimizing csm as a function of cn axis
 */
struct AxisMinimizationFunctor {
  using GroupsType = std::vector<
    std::vector<unsigned>
  >;
  std::reference_wrapper<const PositionCollection> normalizedPositionsRef;
  std::reference_wrapper<const GroupsType> groupsRef;
  unsigned axis_n;

  AxisMinimizationFunctor(
    const PositionCollection& positions,
    const GroupsType& groups,
    const unsigned axisOrder
  ) : normalizedPositionsRef(positions),
      groupsRef(groups),
      axis_n(axisOrder)
  {
    assert(
      temple::all_of(groups,
        [axisOrder](const auto& g) { return g.size() == axisOrder; }
      )
    );
  }

  double evaluate(const Eigen::VectorXd& parameters) {
    const Eigen::Vector3d axis {
      std::sin(parameters(0)) * std::cos(parameters(1)),
      std::sin(parameters(0)) * std::sin(parameters(1)),
      std::cos(parameters(0))
    };

    double csm = 0;
    for(const auto& group : groupsRef.get()) {
      csm += evaluatePermutation(
        normalizedPositionsRef.get(),
        axis,
        axis_n,
        group
      );
    }

    return csm;
  }

  void numericalEvaluation(
    const Eigen::VectorXd& parameters,
    double& value,
    Eigen::Ref<Eigen::VectorXd> gradient,
    Eigen::Ref<Eigen::MatrixXd> hessian
  ) {
    value = evaluate(parameters);
    gradient = temple::optimization::numericalGradient(
      [this](const Eigen::VectorXd& p) -> double {
        return evaluate(p);
      },
      parameters
    );
    hessian = temple::optimization::numericalHessian(
      [this](const Eigen::VectorXd& p) -> double {
        return evaluate(p);
      },
      parameters
    );
  }

  void explicitEvaluation(
    const Eigen::VectorXd& parameters,
    double& value,
    Eigen::Ref<Eigen::VectorXd> gradient,
    Eigen::Ref<Eigen::MatrixXd> hessian
  ) {
    assert(parameters.size() == 2);
    assert(gradient.size() == 2);
    assert(hessian.rows() == 2 && hessian.cols() == 2);

    const double& theta = parameters(0);
    const double& phi = parameters(1);
    const Eigen::Vector3d n {
      std::sin(theta) * std::cos(phi),
      std::sin(theta) * std::sin(phi),
      std::cos(theta)
    };

    assert(std::fabs(n.squaredNorm() - 1) <= 1e-10);

    value = 0;

    Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
    Eigen::Vector3d h = Eigen::Vector3d::Zero();
    const double alpha = 2.0 * M_PI / axis_n;

    const auto& groups = groupsRef.get();
    const auto& normalizedPositions = normalizedPositionsRef.get();

    for(const auto& group : groups) {
      Eigen::Vector3d a, b, c;
      for(unsigned i = 0; i < axis_n; ++i) {
        a.setZero();
        b.setZero();
        c.setZero();

        /* Calculate a, b and c */
        // Terms of j < i
        for(unsigned j = 0; j < i; ++j) {
          const double compoundAngle = (axis_n + j - i) * alpha;
          const double cosine = std::cos(compoundAngle);
          const auto& position = normalizedPositions.col(group.at(j));

          a.noalias() += cosine * position;
          b.noalias() += (1 - cosine) * position;
          c.noalias() += std::sin(compoundAngle) * position;
        }
        // Terms for j == i
        a.noalias() += (1 - static_cast<int>(axis_n)) * normalizedPositions.col(group.at(i));
        // Terms for j > i
        for(unsigned j = i + 1; j < axis_n; ++j) {
          const double compoundAngle = (j - i) * alpha;
          const double cosine = std::cos(compoundAngle);
          const auto& position = normalizedPositions.col(group.at(j));

          a.noalias() += cosine * position;
          b.noalias() += (1 - cosine) * position;
          c.noalias() += std::sin(compoundAngle) * position;
        }

        // Average out all vectors
        a /= axis_n;
        b /= axis_n;
        c /= axis_n;

        // Add value contributions
        value += (a + n.dot(b) * n + n.cross(c)).squaredNorm();

        /* Add contributions to M and h */
        const Eigen::Matrix3d intermediate = a * b.transpose();
        M.noalias() += (
          a * a.transpose()
          - c * c.transpose()
          + intermediate
          + intermediate.transpose()
        );
        h.noalias() += c.cross(a);
      }
    }

    // Calculate gradient components
    const Eigen::Vector3d Mnph = M * n + h;
    const Eigen::Vector3d nPartialTheta {
      std::cos(theta) * std::cos(phi),
      std::cos(theta) * std::sin(phi),
      - std::sin(theta)
    };
    const Eigen::Vector3d nPartialPhi {
      - std::sin(theta) * std::sin(phi),
      std::sin(theta) * std::cos(phi),
      0.0
    };

    gradient(0) = 2 * Mnph.dot(nPartialTheta);
    gradient(1) = 2 * Mnph.dot(nPartialPhi);

    // Calculate hessian components
    const Eigen::Vector3d nPartialThetaPartialPhi {
      - std::cos(theta) * std::sin(phi),
      std::cos(theta) * std::cos(phi),
      0
    };

    // Resolve a 1x1 matrix to its single entry
    auto resolve = [](const Eigen::MatrixXd& matr) -> double {
      assert(matr.rows() == 1 && matr.cols() == 1);
      return matr(0);
    };

    hessian(0, 0) = 2 * (
      resolve(nPartialTheta.transpose() * M * nPartialTheta)
      + Mnph.dot(-n)
    );
    hessian(1, 1) = 2 * (
      resolve(nPartialPhi.transpose() * M * nPartialPhi)
      + Mnph.dot(Eigen::Vector3d {-n(0), -n(1), 0})
    );
    hessian(0, 1) = 2 * (
      resolve(nPartialTheta.transpose() * M * nPartialPhi)
      + Mnph.dot(nPartialThetaPartialPhi)
    );
    hessian(1, 0) = 2 * (
      resolve(nPartialPhi.transpose() * M * nPartialTheta)
      + Mnph.dot(nPartialThetaPartialPhi)
    );
  }

  void operator() (
    const Eigen::VectorXd& parameters,
    double& value,
    Eigen::Ref<Eigen::VectorXd> gradient,
    Eigen::Ref<Eigen::MatrixXd> hessian
  ) {
    numericalEvaluation(parameters, value, gradient, hessian);
    //explicitEvaluation(parameters, value, gradient, hessian);
  }
};

/**
 * @brief Checker for LBFGS minimization
 */
struct AxisMinimizationChecker {
  using FloatType = double;
  using VectorType = Eigen::VectorXd;

  bool shouldContinue(
    const unsigned iteration,
    const FloatType /* value */,
    const VectorType& gradient
  ) {
    return (
      iteration <= 1000
      && gradient.squaredNorm() > 1e-3
    );
  }
};

//! Data struct for axis minimization
struct axis_optimization_t {
  Eigen::Vector3d axis;
  double csm;
};

axis_optimization_t optimize_axis_two_parameters(
  const PositionCollection& normalizedPositions,
  const unsigned n,
  const std::vector<
    std::vector<unsigned>
  >& groups,
  const Eigen::Vector3d& initialAxis
) {
  assert(std::fabs(initialAxis.squaredNorm() - 1) <= 1e-5);

  Eigen::VectorXd parameters (2);
  // Theta from z = cos(theta)
  parameters(0) = std::acos(initialAxis(2));
  // Phi from atan(y / x)
  parameters(1) = std::atan2(initialAxis(1), initialAxis(0));

  AxisMinimizationFunctor functor {
    normalizedPositions,
    groups,
    n
  };

  auto optimizationResult = temple::TrustRegionOptimizer<>::minimize(
    parameters,
    functor,
    AxisMinimizationChecker {}
  );

  if(optimizationResult.iterations >= 1000) {
    throw std::logic_error("Could not minimize axis! Maximum iterations reached.");
  }

  const double theta = parameters(0);
  const double phi = parameters(1);
  Eigen::Vector3d axis {
    std::sin(theta) * std::cos(phi),
    std::sin(theta) * std::sin(phi),
    std::cos(theta)
  };

  return {
    axis,
    100 * optimizationResult.value / (groups.size() * n)
  };
}

struct cn_csm_t {
  double csm;
  std::vector<
    std::vector<unsigned>
  > groups;
};

//! @brief Greedy csm minimization at fixed axis. Not fully permutational.
cn_csm_t fixedAxisGreedy(
  const PositionCollection& normalizedPositions,
  const unsigned n,
  const Eigen::Vector3d& axis
) {
  std::mt19937_64 urbg;
  std::random_device rd;
  urbg.seed(rd());

  assert(n >= 2);
  const unsigned points = normalizedPositions.cols();
  assert(n <= points);

  // Axis vector must be normalized!
  assert(std::fabs(axis.norm() - 1) < 1e-10);

  /* On-axis points should be ignored */
  Eigen::ParametrizedLine<double, 3> axisLine(
    Eigen::Vector3d::Zero(),
    axis
  );

  std::vector<unsigned> offAxisPointIndices;
  offAxisPointIndices.reserve(points);
  for(unsigned i = 0; i < points; ++i) {
    if(axisLine.distance(normalizedPositions.col(i)) > 0.1) {
      offAxisPointIndices.push_back(i);
    }
  }

  const unsigned validPoints = offAxisPointIndices.size();
  /* If there are fewer off-axis points than the rotation order being tested,
   * the element is clearly not present
   */
  if(validPoints < n) {
    cn_csm_t result;
    result.csm = 100;
    return result;
  }

  /* Precalculate the required rotation matrices */
  const double angleRadians = 2 * M_PI / n;
  Eigen::Matrix<double, 3, Eigen::Dynamic> foldMatrices(3, 3 * (n - 1));
  Eigen::Matrix<double, 3, Eigen::Dynamic> unfoldMatrices(3, 3 * (n - 1));
  for(unsigned i = 0; i < n - 1; ++i) {
    foldMatrices.block<3, 3>(0, 3 * i) = Eigen::AngleAxisd((i + 1) * angleRadians, axis).toRotationMatrix();
    unfoldMatrices.block<3, 3>(0, 3 * i) = foldMatrices.block<3, 3>(0, 3 *i).inverse();
  }

  /* Divide off-axis particle indices into groups of size n each */
  const unsigned groups = validPoints / n;
  assert(groups != 0);

  cn_csm_t result;
  result.csm = 100;
  result.groups.resize(groups);

  std::vector<unsigned> groupIndices(validPoints);
  for(unsigned i = 0; i < validPoints; ++i) {
    groupIndices[i] = i / n;
  }

  do {
    /* Groups are subdivided, but unordered. Calculate csm for each subgroup
     * and minimize csm over all permutations
     */
    double groupCollectiveCSM = 0;
    for(unsigned g = 0; g < groups; ++g) {
      // Collect all indices of the current group number
      std::vector<unsigned> group;
      group.reserve(n);
      for(unsigned i = 0; i < validPoints; ++i) {
        if(groupIndices.at(i) == g) {
          group.push_back(offAxisPointIndices.at(i));
        }
      }

      auto calculateCSM = [&](const std::vector<unsigned>& permutation) -> double {
        Eigen::Vector3d averagePoint = normalizedPositions.col(permutation.front());

        /* Fold and average in-place */
        averagePoint.noalias() = normalizedPositions.col(permutation.front());
        for(unsigned i = 1; i < n; ++i) {
          averagePoint.noalias() += foldMatrices.block<3, 3>(0, 3 * (i - 1)) * normalizedPositions.col(permutation.at(i));
        }
        averagePoint /= n;

        /* Calculate CSM while unfolding */
        double csm = (
          normalizedPositions.col(permutation.front())
          - averagePoint
        ).squaredNorm();

        for(unsigned i = 1; i < n; ++i) {
          csm += (
            normalizedPositions.col(permutation.at(i))
            - unfoldMatrices.block<3, 3>(0, 3 * (i - 1)) * averagePoint
          ).squaredNorm();
        }

        return csm;
      };

      assert(group.size() == n);

      std::shuffle(std::begin(group), std::end(group), urbg);

      bool foundBetterAdjacentPermutation = false;
      double currentCSM = calculateCSM(group);
      do {
        foundBetterAdjacentPermutation = false;
        for(unsigned i = 0; i < n && !foundBetterAdjacentPermutation; ++i) {
          for(unsigned j = i + 1; j < n; ++j) {
            std::swap(group.at(i), group.at(j));

            double adjacentCSM = calculateCSM(group);
            if(adjacentCSM < currentCSM) {
              currentCSM = adjacentCSM;
              foundBetterAdjacentPermutation = true;
              break;
            }

            std::swap(group.at(i), group.at(j));
          }
        }
      } while(foundBetterAdjacentPermutation);
      result.groups.at(g) = std::move(group);
      groupCollectiveCSM += currentCSM;
    }

    result.csm = std::min(result.csm, groupCollectiveCSM);
  } while(std::next_permutation(std::begin(groupIndices), std::end(groupIndices)));

  result.csm *= 100.0 / (groups * n);
  return result;
}

//! Minimizes CSM at fixed axis. Fully permutational
double fixedAxis(
  const PositionCollection& normalizedPositions,
  const unsigned n,
  const Eigen::Vector3d& axis
) {
  assert(n >= 2);
  const unsigned points = normalizedPositions.cols();
  assert(n <= points);

  // Axis vector must be normalized!
  assert(std::fabs(axis.norm() - 1) < 1e-10);

  /* On-axis points should be ignored */
  Eigen::ParametrizedLine<double, 3> axisLine(
    Eigen::Vector3d::Zero(),
    axis
  );

  std::vector<unsigned> offAxisPointIndices;
  offAxisPointIndices.reserve(points);
  for(unsigned i = 0; i < points; ++i) {
    if(axisLine.distance(normalizedPositions.col(i)) > 0.1) {
      offAxisPointIndices.push_back(i);
    }
  }

  const unsigned validPoints = offAxisPointIndices.size();
  /* If there are fewer off-axis points than the rotation order being tested,
   * the element is clearly not present
   */
  if(validPoints < n) {
    return 100;
  }

  std::vector<unsigned> groupIndices(validPoints);
  for(unsigned i = 0; i < validPoints; ++i) {
    groupIndices[i] = i / n;
  }

  /* Precalculate the required rotation matrices */
  const double angleRadians = 2 * M_PI / n;
  Eigen::Matrix<double, 3, Eigen::Dynamic> foldMatrices(3, 3 * (n - 1));
  Eigen::Matrix<double, 3, Eigen::Dynamic> unfoldMatrices(3, 3 * (n - 1));
  for(unsigned i = 0; i < n - 1; ++i) {
    foldMatrices.block<3, 3>(0, 3 * i) = Eigen::AngleAxisd((i + 1) * angleRadians, axis).toRotationMatrix();
    unfoldMatrices.block<3, 3>(0, 3 * i) = Eigen::AngleAxisd((i + 1) * angleRadians, -axis).toRotationMatrix();
  }

  /* Divide off-axis particle indices into groups of size n each */
  const unsigned groups = validPoints / n;
  double overallLowestCSM = 100;
  do {
    /* Groups are subdivided, but unordered. Calculate csm for each subgroup
     * and minimize csm there over all permutations
     */
    double groupCollectiveCSM = 0;
    for(unsigned g = 0; g < groups; ++g) {
      // Collect all indices of the current group number
      std::vector<unsigned> group;
      group.reserve(n);
      for(unsigned i = 0; i < validPoints; ++i) {
        if(groupIndices.at(i) == g) {
          group.push_back(offAxisPointIndices.at(i));
        }
      }

      assert(group.size() == n);

      /* All point indices of the current group are ordered ascending right now,
       * so we can iterate through their permutations:
       */
      double lowestCSM = 100;
      do {
        // TODO this can be optimized to avoid foldedPositions and unfoldedPositions completely
        /* Fold */
        Eigen::MatrixXd foldedPositions(3, n);
        foldedPositions.col(0) = normalizedPositions.col(group.front());

        // The first position is unchanged, hence we start with i = 1
        for(unsigned i = 1; i < n; ++i) {
          foldedPositions.col(i) = foldMatrices.block<3, 3>(0, 3 * (i - 1)) * normalizedPositions.col(group.at(i));
        }

        /* Average */
        Eigen::MatrixXd unfoldedPositions(3, n);
        unfoldedPositions.col(0) = foldedPositions.rowwise().sum() / n;

        /* Unfold */
        for(unsigned i = 1; i < n; ++i) {
          unfoldedPositions.col(i) = unfoldMatrices.block<3, 3>(0, 3 * (i - 1)) * unfoldedPositions.col(0);
        }

        /* Calculate CSM */
        double csm = 0;
        for(unsigned i = 0; i < n; ++i) {
          csm += (normalizedPositions.col(group.at(i)) - unfoldedPositions.col(i)).squaredNorm();
        }
        lowestCSM = std::min(lowestCSM, csm);
      } while(std::next_permutation(std::begin(group), std::end(group)));

      groupCollectiveCSM += lowestCSM;
    }

    overallLowestCSM = std::min(overallLowestCSM, groupCollectiveCSM);
  } while(std::next_permutation(std::begin(groupIndices), std::end(groupIndices)));

  return 100 * overallLowestCSM / (groups * n);
}

} // namespace cn

} // namespace csm

InertialMoments principalInertialMoments(
  const PositionCollection& normalizedPositions
) {
  Eigen::Matrix3d inertialMatrix = Eigen::Matrix3d::Zero(3, 3);
  const unsigned N = normalizedPositions.cols();

  for(unsigned i = 0; i < N; ++i) {
    const auto& vec = normalizedPositions.col(i);
    inertialMatrix(0, 0) += vec.y() * vec.y() + vec.z() * vec.z();
    inertialMatrix(1, 1) += vec.x() * vec.x() + vec.z() * vec.z();
    inertialMatrix(2, 2) += vec.x() * vec.x() + vec.y() * vec.y();
    inertialMatrix(1, 0) -= vec.x() * vec.y(); // xy
    inertialMatrix(2, 0) -= vec.x() * vec.z(); // xz
    inertialMatrix(2, 1) -= vec.y() * vec.z(); // yz
  }

  // Decompose the inertial matrix to get principal axes and inertial moments
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> decomposition(inertialMatrix);

  InertialMoments result;
  result.moments = decomposition.eigenvalues();
  result.axes = decomposition.eigenvectors();
  return result;
}

Top standardizeTop(Eigen::Ref<PositionCollection> normalizedPositions) {
  const unsigned N = normalizedPositions.cols();
  assert(N > 1);

  InertialMoments moments = principalInertialMoments(normalizedPositions);

  const unsigned degeneracy = detail::degeneracy(moments.moments);

  /* We can immediately separate out linear cases: If IA << IB = IC and IA ~ 0,
   * then we very likely have a linear molecule on our hands.
   *
   * The principal moments and the corresponding axes are ordered ascending.
   */
  if(moments.moments(0) < 0.1 && degeneracy == 2) {
    return Top::Linear;
  }

  auto rotateEverything = [&](const CoordinateSystem& sourceSystem) {
    const CoordinateSystem defaultCoordinateSystem {};
    assert(sourceSystem.isRightHanded());
    auto R = rotationMatrix(sourceSystem, defaultCoordinateSystem);

    // Rotate coordinates
    for(unsigned i = 0; i < normalizedPositions.cols(); ++i) {
      normalizedPositions.col(i) = R * normalizedPositions.col(i);
    }

    // Rotate inertial moment axes
    for(unsigned i = 0; i < 3; ++i) {
      moments.axes.col(i) = R * moments.axes.col(i);
    }
  };

  if(moments.moments(0) < 0.1 && degeneracy == 2) {
    return Top::Linear;
  }

  if(degeneracy == 1) {
    /* The top is asymmetric. Rotate the axis with the
     * highest moment of inertia to coincide with z, and the one with second most
     * to coincide with x.
     *
     * TODO try to find Cn axes along the inertial moments and then orienting
     * the one with highest order to +z. Fall back to inertial moment ordering
     * if no axes are found.
     */
    CoordinateSystem inertialMomentSystem {
      moments.axes.col(1), // second highest becomes x
      moments.axes.col(2).cross(moments.axes.col(1)) // y = z.cross(x)
    };
    assert(inertialMomentSystem.z.isApprox(moments.axes.col(2), 1e-10));
    rotateEverything(inertialMomentSystem);
    // Make sure rotation went as intended
    assert(moments.axes.col(2).cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));
    assert(moments.axes.col(1).cwiseAbs().isApprox(Eigen::Vector3d::UnitX(), 1e-10));

    return Top::Asymmetric;
  }

  if(degeneracy == 2) {
    /* The top is symmetric. The subsets are:
     * - Oblate (disc): IA = IB < IC
     * - Prolate (rugby football): IA < IB = IC
     *
     * We rotate the unique axis to coincide with z (it's probably the site of
     * the highest-order Cn or Sn, and one of the degenerate axes to coincide
     * with x. There could be a C2 on x.
     *
     * This is most likely rare and should occur only for largely undistorted
     * structures. Perhaps we can flowchart point groups here?
     */
    if(std::fabs((moments.moments(2) - moments.moments(1)) / moments.moments(2)) <= 0.05) {
      // Prolate top. IA is unique
      CoordinateSystem inertialMomentSystem {
        moments.axes.col(1),
        moments.axes.col(2)
      };
      rotateEverything(inertialMomentSystem);
      assert(moments.axes.col(0).cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));
      return Top::Prolate;
    }

    // Oblate top. IC is unique
    CoordinateSystem inertialMomentSystem {
      moments.axes.col(0),
      moments.axes.col(1)
    };
    rotateEverything(inertialMomentSystem);
    assert(moments.axes.col(2).cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));
    return Top::Oblate;
  }

  assert(degeneracy == 3);
  /* The top is spherical (IA = IB = IC).
   *
   * Note that there is no reason to rotate anything on the basis of the axes
   * from the inertial moment analysis since any choice of axes gives the
   * spherical symmetry. We can't use those to rotate the system.
   *
   * Rotate an arbitrary position to +z instead (good for Td
   * and Oh octahedral, less so for Oh cubic and Ih, which should be less
   * common)
   */
  unsigned selectedIndex = 0;
  for(; selectedIndex < N; ++selectedIndex) {
    /* As long as the position isn't close to the centroid and it's not exactly
     * the -z vector, we can rotate it
     */
    if(
      normalizedPositions.col(selectedIndex).norm() > 0.2
      && !normalizedPositions.col(selectedIndex).normalized().isApprox(
        -Eigen::Vector3d::UnitZ(),
        1e-10
      )
    ) {
      break;
    }
  }
  assert(selectedIndex != N);

  // Determine axis of rotation as sum of z and position coordinate
  const Eigen::Vector3d rotationAxis = (
    normalizedPositions.col(selectedIndex).normalized()
    + Eigen::Vector3d::UnitZ()
  ).normalized();
  const Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(M_PI, rotationAxis).toRotationMatrix();

  // Rotate all coordinates
  for(unsigned i = 0; i < N; ++i) {
    normalizedPositions.col(i) = rotationMatrix * normalizedPositions.col(i);
  }
  // Check that everything went as planned
  assert(normalizedPositions.col(selectedIndex).normalized().cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));

  return Top::Spherical;
}

PointGroup flowchart(
  const PositionCollection& normalizedPositions,
  const Top top
) {
  /* Key insights for axis-based flowcharting, if you choose to go that way.
   *
   * An element of symmetry must be parallel or perpendicular to the
   * combination of an adequate number of vectors.
   *
   * Examples:
   * - If you have C3h, and you're testing for a C2 axis perpendicular to the
   *   C3 axis, it will be parallel to a vector i + j (i != j) if n = 6
   *   or to a single i if n = 3.
   * - If you are looking for a mirror plane in C3v, it will be perpendicular
   *   to a vector i - j (i != j)
   * - A C3 axis in a cubic group will be parallel to i if n = 4 or to some
   *   i + j + k (i != j != k) if n > 4
   */

  if(top == Top::Linear) {
    return detail::linear(normalizedPositions);

    /* Do an extra test for collinearity: The rank of the positions matrix is
     * one for linear molecules.
     */
    // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> rankDecomposition(normalizedPositions);
    // rankDecomposition.setThreshold(1e-3);
    // if(rankDecomposition.rank() <= 1) {
    //   // Confirmed linear
    // }
  } else if(top == Top::Spherical) {
    /* Getting a spherical top is rare and will only occur for nearly
     * undistorted structures.  As a consequence, we can flowchart here. Only
     * three point groups have spherical symmetry: Td, Oh, Ih.
     *
     * TODO Maybe a faster way to flowchart is to search for axes along
     * particle positions?
     */
    const double tetrahedralCSM = csm::pointGroup(normalizedPositions, PointGroup::Td).value_or(1000);
    const double octahedralCSM = csm::pointGroup(normalizedPositions, PointGroup::Oh).value_or(10000);
    // TODO icosahedral is also a spherical top!
    return tetrahedralCSM < octahedralCSM ? PointGroup::Td : PointGroup::Oh;
  }

  std::cout << "WARNING: Unimplemented flowchart result\n";
  return PointGroup::C1;
}

  /* Find a main axis and rotate it to z */
  /* Rotate secondary axis to x */
  // for(unsigned axisIndex = 0; axisIndex < 3; ++axisIndex) {
  //   std::cout << "At axis " << axisIndex
  //     << " along " << moments.axes.col(axisIndex).transpose()
  //     << ", inertial moment: " << moments.moments(axisIndex)
  //     << "\n";
  //   for(unsigned n = 2; n <= N; ++n) {
  //     const double axis_csm = csm::cn::fixedAxis(
  //       transformed,
  //       n,
  //       moments.axes.col(axisIndex)
  //     );

  //     const auto csm_greedy = csm::cn::fixedAxisGreedy(
  //       transformed,
  //       n,
  //       moments.axes.col(axisIndex)
  //     );

  //     if(std::fabs(axis_csm - csm_greedy.csm) > 0.1) {
  //       std::cout << "C" << n << " reg and greedy differ strongly: " << axis_csm  << " and " << csm_greedy.csm << "\n";
  //     }

  //     std::cout << "S(C" << n << ") = " << csm_greedy.csm << "\n";

  //     if(axis_csm < 0.1) {
  //       // Cn symmetry csm
  //       const PointGroup cn_point_group = static_cast<PointGroup>(
  //         elements::underlying(PointGroup::C2) + n - 2
  //       );
  //       double symm_csm = csm::pointGroup(
  //         transformed,
  //         cn_point_group
  //       );
  //       std::cout << "C" << n << " point group symmetry csm: " << symm_csm << "\n";
  //     }
  //   }
  // }

} // namespace Symmetry
} // namespace Scine
