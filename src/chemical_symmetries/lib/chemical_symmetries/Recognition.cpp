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
#include "temple/Optimization/SO3NelderMead.h"
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
  assert(p == unfoldMatrices.cols() / 3);

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
  assert(temple::all_of(elementGrouping.groups, [l](const auto& group) { return group.size() == l; }));

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
    /* The source paper proves that the average point is always on some
     * symmetry element and there is no need to do each symmetry element of
     * each group in unfolding
     */
    value += (
      unfoldMatrices.block<3, 3>(0, 3 * elementGrouping.groups.at(i).front()) * averagePoint
      - normalizedPositions.col(particles.at(i))
    ).squaredNorm();
  }

  value *= 100.0 / p;
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
  // TODO remove this stuff when you're happy with the debug
  std::vector<unsigned> bestPermutation;

  do {
    double permutationCSM = calculateCSM(
      normalizedPositions,
      unfoldMatrices,
      foldMatrices,
      particleIndices
    );

    if(permutationCSM < value) {
      value = permutationCSM;
      bestPermutation = particleIndices;
    }
  } while(
    std::next_permutation(
      std::begin(particleIndices),
      std::end(particleIndices)
    )
  );

#ifndef NDEBUG
  static unsigned foldCount = 0;

  std::ofstream outfile("fold" + std::to_string(foldCount) +".xyz");
  const unsigned N = normalizedPositions.cols();
  outfile << (2 * N + 3) << "\n\n";
  outfile << std::fixed << std::setprecision(10);
  /* Regular positions: N */
  for(unsigned i = 0; i < N; ++i) {
    outfile << std::left << std::setw(3) << "H";
    outfile << std::right
      << std::setw(16) << normalizedPositions.col(i).x()
      << std::setw(16) << normalizedPositions.col(i).y()
      << std::setw(16) << normalizedPositions.col(i).z()
      << "\n";
  }
  // Z axis markers: 2
  outfile << std::left << std::setw(3) << "Cl";
  outfile << std::right
    << std::setw(16) << 0.0
    << std::setw(16) << 0.0
    << std::setw(16) << 2.0
    << "\n";
  outfile << std::left << std::setw(3) << "Cl";
  outfile << std::right
    << std::setw(16) << 0.0
    << std::setw(16) << 0.0
    << std::setw(16) << -2.0
    << "\n";
  // Average position: 1
  Eigen::Vector3d averagePoint = Eigen::Vector3d::Zero();
  for(unsigned i = 0; i < N; ++i) {
    averagePoint += foldMatrices.block<3, 3>(0, 3 * i) * normalizedPositions.col(particleIndices.at(i));
  }
  averagePoint /= N;
  outfile << std::left << std::setw(3) << "Br";
  outfile << std::right
      << std::setw(16) << averagePoint.x()
      << std::setw(16) << averagePoint.y()
      << std::setw(16) << averagePoint.z()
      << "\n";
  /* Unfolded positions: N */
  for(unsigned i = 0; i < N; ++i) {
    const Eigen::Vector3d unfolded = unfoldMatrices.block<3, 3>(0, 3 * i) * averagePoint;
    outfile << std::left << std::setw(3) << "F";
    outfile << std::right
        << std::setw(16) << unfolded.x()
        << std::setw(16) << unfolded.y()
        << std::setw(16) << unfolded.z()
        << "\n";
  }

  outfile.close();
  ++foldCount;
#endif

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

  // TODO get rid of bestPermutation

  double value = 1000;
  std::vector<unsigned> bestPermutation;

  do {
    double permutationCSM = calculateCSM(
      normalizedPositions,
      unfoldMatrices,
      foldMatrices,
      particleIndices,
      elementGrouping
    );
    if(permutationCSM < value) {
      value = permutationCSM;
      bestPermutation = particleIndices;
    }
  } while(
    std::next_permutation(
      std::begin(particleIndices),
      std::end(particleIndices)
    )
  );

#ifndef NDEBUG
  static unsigned foldCount = 0;

  std::ofstream outfile("foldgrouped" + std::to_string(foldCount) +".xyz");
  const unsigned N = normalizedPositions.cols();
  const unsigned G = foldMatrices.cols() / 3;
  outfile << (2 * N + 3 + G) << "\n";
  outfile << "CSM  = " << value << "\n";
  outfile << std::fixed << std::setprecision(10);
  /* Regular positions: N */
  // for(unsigned particle : particleIndices) {
  //   outfile << std::left << std::setw(3) << "H";
  //   outfile << std::right
  //     << std::setw(16) << normalizedPositions.col(particle).x()
  //     << std::setw(16) << normalizedPositions.col(particle).y()
  //     << std::setw(16) << normalizedPositions.col(particle).z()
  //     << "\n";
  // }
  for(unsigned i = 0; i < N; ++i) {
    if(temple::makeContainsPredicate(particleIndices)(i)) {
      outfile << std::left << std::setw(3) << "H";
    } else {
      outfile << std::left << std::setw(3) << "C";
    }
    outfile << std::right
      << std::setw(16) << normalizedPositions.col(i).x()
      << std::setw(16) << normalizedPositions.col(i).y()
      << std::setw(16) << normalizedPositions.col(i).z()
      << "\n";
  }
  // Z axis markers: 2
  outfile << std::left << std::setw(3) << "Cl";
  outfile << std::right
    << std::setw(16) << 0.0
    << std::setw(16) << 0.0
    << std::setw(16) << 2.0
    << "\n";
  outfile << std::left << std::setw(3) << "Cl";
  outfile << std::right
    << std::setw(16) << 0.0
    << std::setw(16) << 0.0
    << std::setw(16) << -2.0
    << "\n";
  // Average position: G + 1
  Eigen::Vector3d averagePoint = Eigen::Vector3d::Zero();
  const unsigned p = particleIndices.size();
  const unsigned l = elementGrouping.groups.front().size();
  for(unsigned i = 0; i < p; ++i) {
    const auto& elements = elementGrouping.groups.at(i);
    for(unsigned j = 0; j < l; ++j) {
      const Eigen::Vector3d folded = foldMatrices.block<3, 3>(0, 3 * elements.at(j)) * normalizedPositions.col(bestPermutation.at(i));
      outfile << std::left << std::setw(3) << "N";
      outfile << std::right
        << std::setw(16) << folded.x()
        << std::setw(16) << folded.y()
        << std::setw(16) << folded.z()
        << "\n";
      averagePoint += folded;
    }
  }
  averagePoint /= (p * l);

  outfile << std::left << std::setw(3) << "Br";
  outfile << std::right
      << std::setw(16) << averagePoint.x()
      << std::setw(16) << averagePoint.y()
      << std::setw(16) << averagePoint.z()
      << "\n";
  /* Unfolded positions: N */
  for(unsigned i = 0; i < p; ++i) {
    const auto& elements = elementGrouping.groups.at(i);
    for(unsigned j = 0; j < l; ++j) {
      const Eigen::Vector3d unfolded = unfoldMatrices.block<3, 3>(0, 3 * elements.at(j)) * averagePoint;
      outfile << std::left << std::setw(3) << "F";
      outfile << std::right
          << std::setw(16) << unfolded.x()
          << std::setw(16) << unfolded.y()
          << std::setw(16) << unfolded.z()
          << "\n";
    }
  }

  outfile.close();
  ++foldCount;
#endif

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

  double direct_csm(
    const PositionCollection& positions,
    std::vector<unsigned> particleIndices
  ) const {
    assert(particleIndices.size() == static_cast<decltype(particleIndices.size())>(foldMatrices.cols()) / 3);
    return allSymmetryElements(
      positions,
      unfoldMatrices,
      foldMatrices,
      std::move(particleIndices)
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

    std::cout << "Diophantine constants: " << temple::condense(subdivisionGroupSizes) << " = " << P << "\n";
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
        std::cout << " groups of equal size permutation " << temple::condense(flatGroupMap) << "\n";
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

          std::cout << "  Group of size " << groupSize << " (index " << i << ")\n";

          /* Subpartition csm is minimized over the sub-partitions of
           * groups of equal size
           */
          Partitioner partitioner {multiplier, groupSize};
          double bestSubpartitionCSM = 1000;
          do {
            std::cout << "  Partition: " << temple::condense(partitioner.map()) << "\n";
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

              // TODO move partitionParticles
              const double permutationalGroupCSM = groupedSymmetryElements(
                positions,
                partitionParticles,
                unfoldMatrices,
                foldMatrices,
                npGroup
              );

              std::cout << "   Subpartition: " << temple::condense(partitionParticles) << " CSM = " << permutationalGroupCSM << "\n";

              subpartitionCSM += permutationalGroupCSM;
            }
            subpartitionCSM /= multiplier;
            bestSubpartitionCSM = std::min(bestSubpartitionCSM, subpartitionCSM);
          } while(partitioner.next_partition());

          // Weight by number of points in this subpartition
          partitionCSM += multiplier * groupSize * bestSubpartitionCSM;
        }
        // Average over number of points
        partitionCSM /= P;

        value = std::min(value, partitionCSM);
      } while(std::next_permutation(std::begin(flatGroupMap), std::end(flatGroupMap)));

      std::cout << "Diophantine multipliers: " << temple::condense(subdivisionMultipliers) << ", csm = " << value << "\n";
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

        // TODO reweight?
        const double oneSymmetrizedToOriginCSM = (
          direct_csm(positions, std::move(particleIndices))
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

        // TODO reweight?
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

  double operator() (const Eigen::Matrix3d& rotation) const {
    const PositionCollection rotatedCoordinates = rotation * coordinates;
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

  using MinimizerType = temple::SO3NelderMead<>;
  // Set up the initial simplex to capture asymmetric tops and x/y mixups
  MinimizerType::Parameters simplex;
  simplex[0] = Eigen::Matrix3d::Identity();
  simplex[1] = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitX()).toRotationMatrix();
  simplex[2] = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitY()).toRotationMatrix();
  simplex[3] = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix();

  struct NelderMeadChecker {
    bool shouldContinue(unsigned iteration, double lowestValue, double stddev) const {
      return iteration < 1000 && lowestValue > 1e-3 && stddev > 1e-4;
    }
  };

  auto minimizationResult = MinimizerType::minimize(
    simplex,
    functor,
    NelderMeadChecker {}
  );

#ifndef NDEBUG
  std::cout << "Minimized to " << minimizationResult.value << " in " << minimizationResult.iterations << " iterations.\n";

  static unsigned optimizationNumber = 0;

  auto writeXYZ = [](const std::string& filename, const PositionCollection& positions) {
    std::ofstream outfile(filename);
    const unsigned N = positions.cols();
    outfile << (N + 2) << "\n\n";
    outfile << std::fixed << std::setprecision(10);
    for(unsigned i = 0; i < N; ++i) {
      outfile << std::left << std::setw(3) << "H";
      outfile << std::right
        << std::setw(16) << positions.col(i).x()
        << std::setw(16) << positions.col(i).y()
        << std::setw(16) << positions.col(i).z()
        << "\n";
    }
    outfile << std::left << std::setw(3) << "F";
    outfile << std::right
      << std::setw(16) << 0.0
      << std::setw(16) << 0.0
      << std::setw(16) << 2.0
      << "\n";
    outfile << std::left << std::setw(3) << "F";
    outfile << std::right
      << std::setw(16) << 0.0
      << std::setw(16) << 0.0
      << std::setw(16) << -2.0
      << "\n";
    outfile.close();
  };

  writeXYZ(
    std::to_string(optimizationNumber)+".xyz",
    simplex.at(minimizationResult.minimalIndex) * normalizedPositions
  );

  ++optimizationNumber;
#endif

  return minimizationResult.value;
}

double element(
  const PositionCollection& normalizedPositions,
  elements::Rotation rotation
) {
  assert(std::fabs(rotation.axis.norm() - 1) < 1e-10);
  assert(rotation.n >= 2);
  const unsigned P = normalizedPositions.cols();

  if(rotation.n > P) {
    return 100;
  }

  Eigen::ParametrizedLine<double, 3> axisLine(
    Eigen::Vector3d::Zero(),
    rotation.axis
  );

  std::vector<unsigned> diophantineConstants {rotation.n, 1};
  if(!diophantine::has_solution(diophantineConstants, P)) {
    return 100;
  }

  std::vector<unsigned> diophantineMultipliers;
  if(!diophantine::first_solution(diophantineMultipliers, diophantineConstants, P)) {
    throw std::logic_error("Diophantine failure! Couldn't find first solution");
  }

  /* Precalculate fold and unfold matrices */
  Eigen::Matrix<double, 3, Eigen::Dynamic> foldMatrices(3, 3 * (rotation.n - 1));
  Eigen::Matrix<double, 3, Eigen::Dynamic> unfoldMatrices(3, 3 * (rotation.n - 1));
  for(unsigned i = 0; i < rotation.n - 1; ++i) {
    rotation.power = (i + 1);
    foldMatrices.block<3, 3>(0, 3 * i) = rotation.matrix();
    unfoldMatrices.block<3, 3>(0, 3 * i) = foldMatrices.block<3, 3>(0, 3 * i).inverse();
  }

  auto calculateBestPermutationCSM = [&](std::vector<unsigned> particlePermutation) -> double {
    assert(std::is_sorted(std::begin(particlePermutation), std::end(particlePermutation)));
    const unsigned p = particlePermutation.size();

    double minimalPermutationCSM = 1000;
    do {
      /* Fold and average */
      Eigen::Vector3d averagePoint = normalizedPositions.col(particlePermutation.front());
      for(unsigned i = 1; i < p; ++i) {
        averagePoint += foldMatrices.block<3, 3>(0, 3 * (i - 1)) * normalizedPositions.col(particlePermutation.at(i));
      }
      averagePoint /= p;

      /* Unfold and calculate CSM */
      double csm = (normalizedPositions.col(particlePermutation.front()) - averagePoint).squaredNorm();
      for(unsigned i = 1; i < p; ++i) {
        csm += (
          unfoldMatrices.block<3, 3>(0, 3 * (i - 1)) * averagePoint
          - normalizedPositions.col(particlePermutation.at(i))
        ).squaredNorm();
      }
      csm /= p;

      minimalPermutationCSM = std::min(minimalPermutationCSM, csm);
    } while(std::next_permutation(std::begin(particlePermutation), std::end(particlePermutation)));

    return minimalPermutationCSM;
  };

  double value = 1000;
  do {
    // Handle case that all points are symmetrized
    if(diophantineMultipliers.front() == 0) {
      double allAxisSymmetrizedCSM = 0;
      for(unsigned i = 0; i < P; ++i) {
        allAxisSymmetrizedCSM += axisLine.squaredDistance(normalizedPositions.col(i));
      }
      allAxisSymmetrizedCSM /= P;
      value = std::min(value, allAxisSymmetrizedCSM);
      continue;
    }

    std::vector<unsigned> partitionOrAxisSymmetrize;
    partitionOrAxisSymmetrize.reserve(P);
    partitionOrAxisSymmetrize.resize(rotation.n * diophantineMultipliers.front(), 0);
    partitionOrAxisSymmetrize.resize(P, 1);

    double diophantineCSM = 1000;
    do {
      double permutationCSM = 0;
      /* Collect indices to partition */
      std::vector<unsigned> indicesToPartition;
      for(unsigned i = 0; i < P; ++i) {
        if(partitionOrAxisSymmetrize.at(i) == 0) {
          indicesToPartition.push_back(i);
        }
      }
      assert(!indicesToPartition.empty());

      /* Perform partitioning */
      Partitioner partitioner {diophantineMultipliers.front(), rotation.n};
      assert(diophantineMultipliers.front() * rotation.n == indicesToPartition.size());
      double bestPartitionCSM = 1000;
      do {
        double partitionCSM = 0;
        for(auto&& partitionIndices : partitioner.partitions()) {
          auto partitionParticleIndices = temple::map(partitionIndices, temple::functor::at(indicesToPartition));
          const double subpartitionCSM = calculateBestPermutationCSM(partitionParticleIndices);
          partitionCSM += subpartitionCSM;
        }
        partitionCSM /= diophantineMultipliers.front();
        bestPartitionCSM = std::min(bestPartitionCSM, partitionCSM);
      } while(partitioner.next_partition());
      permutationCSM += rotation.n * diophantineMultipliers.front() * bestPartitionCSM;

      /* Add indices to axis symmetrize contributions to CSM */
      for(unsigned i = 0; i < P; ++i) {
        if(partitionOrAxisSymmetrize.at(i) == 1) {
          permutationCSM += axisLine.squaredDistance(normalizedPositions.col(i));
        }
      }

      permutationCSM /= P;
      diophantineCSM = std::min(diophantineCSM, permutationCSM);
    } while(std::next_permutation(std::begin(partitionOrAxisSymmetrize), std::end(partitionOrAxisSymmetrize)));
    value = std::min(value, diophantineCSM);
  } while(diophantine::next_solution(diophantineMultipliers, diophantineConstants, P));
  return 100 * value;
}

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
    // The top is linear: If IA << IB = IC and IA ~ 0. We rotate IA to z
    const CoordinateSystem inertialMomentSystem {
      moments.axes.col(1),
      moments.axes.col(2)
    };
    rotateEverything(inertialMomentSystem);
    assert(moments.axes.col(0).cwiseAbs().isApprox(Eigen::Vector3d::UnitZ(), 1e-10));
    return Top::Linear;
  }

  if(degeneracy == 1) {
    /* The top is asymmetric. Rotate the axis with the
     * highest moment of inertia to coincide with z, and the one with second most
     * to coincide with x.
     *
     * To better define orientation, we could look for Cn axes. This is done in
     * another function. No need to burden this function with that here.
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
    // Calculate Ray's asymmetry parameter
    const double A = 1 / moments.moments(0);
    const double B = 1 / moments.moments(1);
    const double C = 1 / moments.moments(2);
    const double kappa = (2 * B - A - C) / (A - C);
    assert(-1 <= kappa && kappa <= 1);
    if(kappa < 0) {
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

unsigned reorientAsymmetricTop(Eigen::Ref<PositionCollection> normalizedPositions) {
  const unsigned P = normalizedPositions.cols();
  const auto& axes = Eigen::Matrix3d::Identity();

  struct AxisBest {
    unsigned order = 1;
    double csm = 1; // This functions much like a detection threshold below
    unsigned axisIndex;

    AxisBest(unsigned index) : axisIndex(index) {}

    bool operator < (const AxisBest& other) const {
      return order > other.order;
    }
  };

  auto orderedAxisBest = temple::sort(
    temple::map(
      temple::iota<unsigned>(3),
      [&](const unsigned axisIndex) -> AxisBest {
        const Eigen::Vector3d axis = axes.col(axisIndex);
        AxisBest best {axisIndex};
        for(unsigned n = 2; n <= P; ++n) {
          const double axisCSM = csm::element(normalizedPositions, elements::Rotation::Cn(axis, n));
          if(axisCSM < best.csm) {
            best.order = n;
            best.csm = axisCSM;
          }
        }
        return best;
      }
    )
  );

  if(orderedAxisBest.front().order > 1) {
    /* Only mess with the coordinate frame if any sort of axis was found.
     * We want the second-highest order axis on x, highest order axis along z,
     * doesn't really matter if +z or -z
     */
    const CoordinateSystem highestOrderSystem {
      axes.col(orderedAxisBest.at(1).axisIndex),
      axes.col(orderedAxisBest.back().axisIndex)
    };

    normalizedPositions = rotationMatrix(highestOrderSystem, {}) * normalizedPositions;
  }

  return orderedAxisBest.front().order;
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
    const double octahedralCSM = csm::pointGroup(normalizedPositions, PointGroup::Oh).value_or(1000);
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
