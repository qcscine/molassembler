#include "shapes/ContinuousMeasures.h"

#include "shapes/Diophantine.h"
#include "shapes/Partitioner.h"
#include "shapes/Symmetries.h"

#include "temple/Adaptors/Iota.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"
#include "temple/Optimization/SO3NelderMead.h"
#include "temple/constexpr/Numeric.h"

#include "boost/math/tools/minima.hpp"
#include <Eigen/Eigenvalues>

namespace Scine {
namespace Symmetry {
namespace continuous {

using Matrix = Eigen::Matrix<double, 3, Eigen::Dynamic>;

template<typename Derived>
bool centroidIsZero(const Eigen::MatrixBase<Derived>& a) {
  assert(a.rows() == 3);
  return (a.rowwise().sum() / a.cols()).squaredNorm() < 1e-8;
}

template<typename DerivedA, typename DerivedB>
Eigen::Matrix3d fitQuaternion(const Eigen::MatrixBase<DerivedA>& stator, const Eigen::MatrixBase<DerivedB>& rotor) {
  assert(centroidIsZero(stator));
  assert(centroidIsZero(rotor));

  Eigen::Matrix4d b = Eigen::Matrix4d::Zero();
  // generate decomposable matrix per atom and add them
  for (int i = 0; i < rotor.cols(); i++) {
    auto& rotorCol = rotor.col(i);
    auto& statorCol = stator.col(i);

    Eigen::Matrix4d a = Eigen::Matrix4d::Zero();
    a.block<1, 3>(0, 1) = (rotorCol - statorCol).transpose();
    a.block<3, 1>(1, 0) = statorCol - rotorCol;
    a.block<3, 3>(1, 1) = Eigen::Matrix3d::Identity().rowwise().cross(statorCol + rotorCol);
    b += a.transpose() * a;
  }

  // Decompose b
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigensolver(b);

  // Do not allow improper rotation
  const Eigen::Vector4d& q = eigensolver.eigenvectors().col(0);
  return Eigen::Quaterniond(q[0], q[1], q[2], q[3]).toRotationMatrix();
}

Eigen::Matrix3d fitQuaternion(
  const PositionCollection& stator,
  const PositionCollection& rotor,
  const std::unordered_map<unsigned, unsigned>& p
) {
  assert(centroidIsZero(stator));
  assert(centroidIsZero(rotor));

  Eigen::Matrix4d b = Eigen::Matrix4d::Zero();
  // generate decomposable matrix per atom and add them
  for(auto& iterPair : p) {
    auto& statorCol = stator.col(iterPair.first);
    auto& rotorCol = rotor.col(iterPair.second);

    Eigen::Matrix4d a = Eigen::Matrix4d::Zero();
    a.block<1, 3>(0, 1) = (rotorCol - statorCol).transpose();
    a.block<3, 1>(1, 0) = statorCol - rotorCol;
    a.block<3, 3>(1, 1) = Eigen::Matrix3d::Identity().rowwise().cross(statorCol + rotorCol);
    b += a.transpose() * a;
  }

  // Decompose b
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigensolver(b);

  // Do not allow improper rotation
  const Eigen::Vector4d& q = eigensolver.eigenvectors().col(0);
  return Eigen::Quaterniond(q[0], q[1], q[2], q[3]).toRotationMatrix();
}

bool isNormalized(const PositionCollection& positions) {
  /* All vectors are shorter or equally long as a unit vector
   * and at least one vector is as long as a unit vector
   */
  const unsigned P = positions.cols();
  return (
    temple::all_of(
      temple::iota<unsigned>(P),
      [&](const unsigned i) -> bool {
        return positions.col(i).norm() <= (1 + 1e-5);
      }
    ) && temple::any_of(
      temple::iota<unsigned>(P),
      [&](const unsigned i) -> bool {
        return std::fabs(positions.col(i).norm() - 1) < 1e-5;
      }
    )
  );
}

PositionCollection normalize(const PositionCollection& positions) {
  // Translate the origin to the average position
  const Eigen::Vector3d center = positions.rowwise().sum() / positions.cols();
  PositionCollection transformed = positions.colwise() - center;

  // Rescale all distances so that the longest is a unit vector
  transformed /= std::sqrt(transformed.colwise().squaredNorm().maxCoeff());
  assert(isNormalized(transformed));
  return transformed;
}

namespace fixed {

double element(
  const PositionCollection& normalizedPositions,
  const elements::Rotation& rotation
) {
  assert(rotation.power == 1);
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
  auto cumulativeRotation = rotation;
  for(unsigned i = 0; i < rotation.n - 1; ++i) {
    foldMatrices.block<3, 3>(0, 3 * i) = cumulativeRotation.matrix();
    unfoldMatrices.block<3, 3>(0, 3 * i) = foldMatrices.block<3, 3>(0, 3 * i).inverse();

    ++cumulativeRotation.power;
    cumulativeRotation.reflect xor_eq rotation.reflect;
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
      /* Collect indices to partition, axis symmetrize the rest */
      std::vector<unsigned> indicesToPartition;
      for(unsigned i = 0; i < P; ++i) {
        if(partitionOrAxisSymmetrize.at(i) == 0) {
          indicesToPartition.push_back(i);
        } else {
          permutationCSM += axisLine.squaredDistance(normalizedPositions.col(i));
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
      permutationCSM /= P;
      diophantineCSM = std::min(diophantineCSM, permutationCSM);
    } while(std::next_permutation(std::begin(partitionOrAxisSymmetrize), std::end(partitionOrAxisSymmetrize)));
    value = std::min(value, diophantineCSM);
  } while(diophantine::next_solution(diophantineMultipliers, diophantineConstants, P));

  return 100 * value;
}

double element(
  const PositionCollection& normalizedPositions,
  const elements::Reflection& reflection
) {
  /* TODO
   * - calculateReflectionCSM could be memoized across its a,b arguments
   *   also try memoizing other elements' basic calculation fns
   */

  const unsigned P = normalizedPositions.cols();
  assert(P > 0);

  Eigen::Hyperplane<double, 3> plane(reflection.normal, 0);
  const Eigen::Matrix3d reflectMatrix = reflection.matrix();
  assert(reflectMatrix.inverse().isApprox(reflectMatrix, 1e-10));

  const std::vector<unsigned> diophantineConstants {2, 1};
  assert(diophantine::has_solution(diophantineConstants, P));
  std::vector<unsigned> diophantineMultipliers;
  if(!diophantine::first_solution(diophantineMultipliers, diophantineConstants, P)) {
    throw std::logic_error("Diophantine failure! Couldn't find first solution");
  }

  auto calculateReflectionCSM = [](
    const unsigned a,
    const unsigned b,
    const Eigen::Matrix3d& R,
    const PositionCollection& positions
  ) -> double {
    /* Average */
    const Eigen::Vector3d average = (positions.col(a) + R * positions.col(b)) / 2;
    return (
      (positions.col(a) - average).squaredNorm()
      + ((R * average) - positions.col(b)).squaredNorm()
    ) / 2;
  };

  double value = 1000;
  do {
    if(diophantineMultipliers.front() == 0) {
      /* All points are symmetrized to the plane */
      const double csm = temple::sum(
        temple::adaptors::transform(
          temple::adaptors::range(P),
          [&](const unsigned i) -> double {
            return std::pow(plane.absDistance(normalizedPositions.col(i)), 2);
          }
        )
      ) / P;
      value = std::min(value, csm);
      continue;
    }

    std::vector<unsigned> pairOrPlaneSymmetrize;
    pairOrPlaneSymmetrize.reserve(P);
    pairOrPlaneSymmetrize.resize(2 * diophantineMultipliers.front(), 0);
    pairOrPlaneSymmetrize.resize(P, 1);

    double diophantineCSM = 1000;
    do {
      double permutationCSM = 0;
      /* Collect indices to partition */
      std::vector<unsigned> indicesToPartition;
      for(unsigned i = 0; i < P; ++i) {
        if(pairOrPlaneSymmetrize.at(i) == 0) {
          indicesToPartition.push_back(i);
        } else {
          permutationCSM += std::pow(plane.absDistance(normalizedPositions.col(i)), 2);
        }
      }

      /* Perform partitioning into groups of size two */
      Partitioner partitioner {diophantineMultipliers.front(), 2};
      assert(diophantineMultipliers.front() * 2 == indicesToPartition.size());
      double bestPartitionCSM = 1000;
      do {
        double partitionCSM = 0;
        for(auto&& partitionIndices : partitioner.partitions()) {
          assert(partitionIndices.size() == 2);
          const unsigned i = indicesToPartition.at(partitionIndices.front());
          const unsigned j = indicesToPartition.at(partitionIndices.back());
          const double subpartitionCSM = calculateReflectionCSM(i, j, reflectMatrix, normalizedPositions);
          partitionCSM += subpartitionCSM;
        }
        partitionCSM /= 2;
        bestPartitionCSM = std::min(bestPartitionCSM, partitionCSM);
      } while(partitioner.next_partition());
      permutationCSM += 2 * diophantineMultipliers.front() * bestPartitionCSM;
      permutationCSM /= P;
      diophantineCSM = std::min(diophantineCSM, permutationCSM);
    } while(std::next_permutation(std::begin(pairOrPlaneSymmetrize), std::end(pairOrPlaneSymmetrize)));
    value = std::min(value, diophantineCSM);
  } while(diophantine::next_solution(diophantineMultipliers, diophantineConstants, P));

  return 100 * value;
}

double Cinf(
  const PositionCollection& normalizedPositions,
  const Eigen::Vector3d& axis
) {
  Eigen::ParametrizedLine<double, 3> axisLine(
    Eigen::Vector3d::Zero(),
    axis
  );
  const unsigned P = normalizedPositions.cols();
  double value = 0;
  for(unsigned i = 0; i < P; ++i) {
    value += axisLine.squaredDistance(normalizedPositions.col(i));
  }
  return 100 * value / P;
}

} // namespace fixed

std::pair<double, elements::Rotation> element(
  const PositionCollection& normalizedPositions,
  elements::Rotation rotation
) {
  using MinimizerType = temple::SO3NelderMead<>;

  MinimizerType::Parameters simplex;
  simplex.at(0) = Eigen::Matrix3d::Identity();
  simplex.at(1) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitX()).toRotationMatrix();
  simplex.at(2) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitY()).toRotationMatrix();
  simplex.at(3) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix();

  struct NelderMeadChecker {
    bool shouldContinue(unsigned iteration, double lowestValue, double stddev) const {
      return iteration < 1000 && lowestValue > 1e-3 && stddev > 1e-4;
    }
  };

  auto minimizationResult = MinimizerType::minimize(
    simplex,
    [&](const Eigen::Matrix3d& R) -> double {
      return fixed::element(R * normalizedPositions, rotation);
    },
    NelderMeadChecker {}
  );

  rotation.axis = simplex.at(minimizationResult.minimalIndex).inverse() * rotation.axis;

  return {
    minimizationResult.value,
    rotation
  };
}


std::pair<double, elements::Reflection> element(
  const PositionCollection& normalizedPositions,
  elements::Reflection reflection
) {
  using MinimizerType = temple::SO3NelderMead<>;

  MinimizerType::Parameters simplex;
  simplex.at(0) = Eigen::Matrix3d::Identity();
  simplex.at(1) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitX()).toRotationMatrix();
  simplex.at(2) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitY()).toRotationMatrix();
  simplex.at(3) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix();

  struct NelderMeadChecker {
    bool shouldContinue(unsigned iteration, double lowestValue, double stddev) const {
      return iteration < 1000 && lowestValue > 1e-3 && stddev > 1e-4;
    }
  };

  auto minimizationResult = MinimizerType::minimize(
    simplex,
    [&](const Eigen::Matrix3d& R) -> double {
      return fixed::element(R * normalizedPositions, reflection);
    },
    NelderMeadChecker {}
  );

  return {
    minimizationResult.value,
    elements::Reflection {
      simplex.at(minimizationResult.minimalIndex).inverse() * reflection.normal
    }
  };
}

double element(
  const PositionCollection& normalizedPositions,
  const elements::Inversion& /* inversion */
) {
  auto calculateInversionCSM = [&](
    const unsigned a,
    const unsigned b,
    const PositionCollection& positions
  ) -> double {
    const Eigen::Vector3d average = (positions.col(a) - positions.col(b)) / 2;
    return (
      (positions.col(a) - average).squaredNorm()
      + (positions.col(b) + average).squaredNorm()
    );
  };

  // If the number of points is even, we can go fast
  const unsigned P = normalizedPositions.cols();
  if(P % 2 == 0) {
    double bestPartitionCSM = 1000;
    Partitioner partitioner {P / 2, 2};
    do {
      double partitionCSM = 0;
      for(auto&& partitionIndices : partitioner.partitions()) {
        assert(partitionIndices.size() == 2);
        const unsigned i = partitionIndices.front();
        const unsigned j = partitionIndices.back();
        partitionCSM += calculateInversionCSM(i, j, normalizedPositions);
      }
      bestPartitionCSM = std::min(bestPartitionCSM, partitionCSM);
    } while(partitioner.next_partition());
    return 100 * bestPartitionCSM / P;
  }

  /* If the number of points is not even, then one point may be symmetrized to
   * the origin, but all others must be paired off and CSM calculated like for
   * Reflection.
   */
  double bestCSM = 1000;

  for(unsigned excludedIndex = 0; excludedIndex < P; ++excludedIndex) {
    // i is the index to exclude from partitioning
    std::vector<unsigned> indicesToPartition;
    indicesToPartition.reserve(P - 1);
    for(unsigned j = 0; j < excludedIndex; ++j) {
      indicesToPartition.push_back(j);
    }
    for(unsigned j = excludedIndex + 1; j < P; ++j) {
      indicesToPartition.push_back(j);
    }

    double bestPartitionCSM = 1000;
    Partitioner partitioner {P / 2, 2};
    do {
      double partitionCSM = 0;
      for(auto&& partitionIndices : partitioner.partitions()) {
        assert(partitionIndices.size() == 2);
        const unsigned i = indicesToPartition.at(partitionIndices.front());
        const unsigned j = indicesToPartition.at(partitionIndices.back());
        partitionCSM += calculateInversionCSM(i, j, normalizedPositions);
      }
      bestPartitionCSM = std::min(bestPartitionCSM, partitionCSM);
    } while(partitioner.next_partition());

    const double csmValue = bestPartitionCSM + normalizedPositions.col(excludedIndex).squaredNorm();
    bestCSM = std::min(bestCSM, csmValue);
  }

  return 100 * bestCSM / P;
}

double Cinf(const PositionCollection& normalizedPositions) {
  struct Functor {
    const PositionCollection& coordinates;
    Functor(const PositionCollection& normalizedPositions) : coordinates(normalizedPositions) {}
    double operator() (const Eigen::Matrix3d& rotation) const {
      const PositionCollection rotatedCoordinates = rotation * coordinates;
      return fixed::Cinf(rotatedCoordinates, Eigen::Vector3d::UnitZ());
    }
  };

  Functor functor {normalizedPositions};

  using MinimizerType = temple::SO3NelderMead<>;

  // Set up the initial simplex to capture asymmetric tops and x/y mixups
  MinimizerType::Parameters simplex;
  simplex.at(0) = Eigen::Matrix3d::Identity();
  simplex.at(1) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitX()).toRotationMatrix();
  simplex.at(2) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitY()).toRotationMatrix();
  simplex.at(3) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix();

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

  return minimizationResult.value;
}

double calculateCSM(
  const PositionCollection& normalizedPositions,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  const std::vector<unsigned>& particles
) {
  const unsigned p = particles.size();
  assert(p == unfoldMatrices.cols() / 3);

  /* Fold points and average */
  Eigen::Vector3d averagePoint = Eigen::Vector3d::Zero();
  for(unsigned i = 0; i < p; ++i) {
    averagePoint += foldMatrices.block<3, 3>(0, 3 * i) * normalizedPositions.col(particles.at(i));
  }
  averagePoint /= p;

  /* Unfold points */
  double value = 0;
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
  double value = 0;
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


/*! @brief Minimizes CSM for a point group, case: G = P
 *
 * This minimizes the continuous symmetry measure for the case that the number
 * of group symmetry elements matches the number of particles.
 */
double allSymmetryElements(
  const PositionCollection& normalizedPositions,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  std::vector<unsigned> particleIndices
) {
  assert(particleIndices.size() == static_cast<std::size_t>(unfoldMatrices.cols()) / 3);

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

/*! @brief Minimizes CSM for a point group, case G = l * P
 *
 * This minimizes the continuous symmetry measure for the case that the number
 * of group symmetry elements is a multiple l of the number of particles.
 *
 * @todo consider particleIndices by-ref
 */
double groupedSymmetryElements(
  const PositionCollection& normalizedPositions,
  std::vector<unsigned> particleIndices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& unfoldMatrices,
  const Eigen::Matrix<double, 3, Eigen::Dynamic>& foldMatrices,
  const std::vector<elements::ElementGrouping>& elementGroupings
) {
  /* The number of groups in element grouping must match the number of particles
   * permutated here
   */
  assert(std::is_sorted(std::begin(particleIndices), std::end(particleIndices)));

  double value = 1000;
  do {
    double minPermutation = 1000;
    for(const auto& grouping : elementGroupings) {
      minPermutation = std::min(
        minPermutation,
        calculateCSM(
          normalizedPositions,
          unfoldMatrices,
          foldMatrices,
          particleIndices,
          grouping
        )
      );
    }
    value = std::min(value, minPermutation);
  } while(
    std::next_permutation(
      std::begin(particleIndices),
      std::end(particleIndices)
    )
  );

  return value;
}

struct OrientationCSMFunctor {
  using MatrixType = Eigen::Matrix<double, 3, Eigen::Dynamic>;

  const PositionCollection& coordinates;
  const MatrixType unfoldMatrices;
  const MatrixType foldMatrices;
  const elements::NPGroupingsMapType npGroups;

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

  double diophantine_csm(
    const PositionCollection& positions,
    const std::vector<unsigned>& subdivisionGroupSizes,
    const std::vector<unsigned>& particleIndices
  ) const {
    const unsigned P = particleIndices.size();
    const unsigned G = foldMatrices.cols() / 3;
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
          const unsigned multiplier = subdivisionMultipliers.at(i);
          // Skip groups with multiplier zero
          if(multiplier == 0) {
            continue;
          }

          const unsigned groupSize = subdivisionGroupSizes.at(i);
          const auto& sameSizeParticleIndices = sameSizeIndexGroups.at(i);

          // Catch group sizes equal to the number of symmetry elements
          if(groupSize == G) {
            Partitioner partitioner {multiplier, G};
            double bestPartitionCSM = 1000;
            do {
              double partitionCSMSum = 0;
              for(auto&& partitionIndices : partitioner.partitions()) {
                auto partitionParticles = temple::map(
                  partitionIndices,
                  [&](const unsigned indexOfParticle) -> unsigned {
                    return particleIndices.at(
                      sameSizeParticleIndices.at(indexOfParticle)
                    );
                  }
                );

                const double subpartitionCSM = allSymmetryElements(
                  positions,
                  unfoldMatrices,
                  foldMatrices,
                  partitionParticles
                );

                partitionCSMSum += subpartitionCSM;
              }
              partitionCSMSum /= multiplier;
              bestPartitionCSM = std::min(bestPartitionCSM, partitionCSMSum);
            } while(partitioner.next_partition());

            partitionCSM += multiplier * groupSize * bestPartitionCSM;
            continue;
          }

          const auto& npGroup = npGroups.at(groupSize);


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

              // TODO move partitionParticles
              const double permutationalGroupCSM = groupedSymmetryElements(
                positions,
                partitionParticles,
                unfoldMatrices,
                foldMatrices,
                npGroup
              );

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
    } while(diophantine::next_solution(subdivisionMultipliers, subdivisionGroupSizes, P));

    return value;
  }

  double csm(const PositionCollection& positions) const {
    const unsigned P = positions.cols();
    const unsigned G = foldMatrices.cols() / 3;

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

    throw std::logic_error("You shouldn't even instantiate this type if you know that you cannot calculate a CSM");
  }

  double operator() (const Eigen::Matrix3d& rotation) const {
    const PositionCollection rotatedCoordinates = rotation * coordinates;
    return csm(rotatedCoordinates);
  }
};

double pointGroup(
  const PositionCollection& normalizedPositions,
  const PointGroup group
) {
  assert(isNormalized(normalizedPositions));

  // Special-case Cinfv
  if(group == PointGroup::Cinfv) {
    return Cinf(normalizedPositions);
  }

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

  if(!diophantine::has_solution(subdivisionGroupSizes, P)) {
    throw std::logic_error("Cannot calculate a CSM for this number of points and this point group. This is most likely an implementation error.");
  }

  OrientationCSMFunctor functor {normalizedPositions, group};

  using MinimizerType = temple::SO3NelderMead<>;
  // Set up the initial simplex to capture asymmetric tops and x/y mixups
  MinimizerType::Parameters simplex;
  simplex.at(0) = Eigen::Matrix3d::Identity();
  simplex.at(1) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitX()).toRotationMatrix();
  simplex.at(2) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitY()).toRotationMatrix();
  simplex.at(3) = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix();

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

  return minimizationResult.value;
}

double shapeFaithfulPaperImplementation(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
  assert(isNormalized(normalizedPositions));
  const unsigned N = normalizedPositions.cols();

  if(N != size(shape) + 1) {
    throw std::logic_error("Mismatched number of positions between supplied coordinates and shape!");
  }

  auto permutation = temple::iota<unsigned>(N);

  // Add the origin
  Matrix shapeCoordinates(3, N);
  shapeCoordinates.block(0, 0, 3, N - 1) = symmetryData().at(shape).coordinates;
  shapeCoordinates.col(N - 1) = Eigen::Vector3d::Zero();
  // Normalize the coordinates
  shapeCoordinates = normalize(shapeCoordinates);

  constexpr double scalingLowerBound = 0.5;
  constexpr double scalingUpperBound = 1.1;

  double permutationalMinimum = std::numeric_limits<double>::max();

  Eigen::Matrix<double, 3, Eigen::Dynamic> permutedShape(3, N);
  std::vector<unsigned> bestPermutation;
  Eigen::Matrix3d bestRotationMatrix;
  do {
    // Construct a permuted shape positions matrix
    for(unsigned i = 0; i < N; ++i) {
      permutedShape.col(permutation.at(i)) = shapeCoordinates.col(i);
    }

    // Perform a quaternion fit to minimize square norm difference over rotation
    auto rotationMatrix = fitQuaternion(normalizedPositions, permutedShape);
    permutedShape = rotationMatrix * permutedShape;

    // Minimize over isotropic scaling factor
    auto scalingMinimizationResult = boost::math::tools::brent_find_minima(
      [&](const double scaling) -> double {
        return (normalizedPositions - scaling * permutedShape).colwise().squaredNorm().sum();
      },
      scalingLowerBound,
      scalingUpperBound,
      std::numeric_limits<double>::digits
    );

    permutationalMinimum = std::min(permutationalMinimum, scalingMinimizationResult.second);
  } while(std::next_permutation(std::begin(permutation), std::end(permutation)));

  const double normalization = normalizedPositions.colwise().squaredNorm().sum();
  return  100 * permutationalMinimum / normalization;
}

double shapeAlternateImplementation(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
  /* There is a bit of a conundrum: The paper says
   * - For each index mapping
   *   - Minimize over rotation
   *   - Minimize over scaling
   *
   * But it's a lot faster to
   * - Minimize over rotation (while minimizing CShapeM over all permutations)
   * - Minimize over scaling (using best permutation from rotation step)
   *
   * And to me it's not immediately apparent why this should be worse, especially
   * considering that minimizing over permutations while calculating CShapeM
   * should be smooth. Maybe it has local minima that the paper procedure wouldn't?
   *
   * But it's odd. The paper suggests pre-pairing off vertices to reduce cost "in
   * most cases". Not sure which is more dangerous (pre-pairing prior to any
   * minimizations or reversing minimization order and reusing pairing).
   *
   * So here we speed up the faithful implementation by minimizing over rotation,
   * remembering the best rotation matrix and minimizing over scaling outside
   * of the permutational loop
   */

  assert(isNormalized(normalizedPositions));
  const unsigned N = normalizedPositions.cols();

  if(N != size(shape) + 1) {
    throw std::logic_error("Mismatched number of positions between supplied coordinates and shape!");
  }

  auto permutation = temple::iota<unsigned>(N);

  // Add the origin
  Matrix shapeCoordinates(3, N);
  shapeCoordinates.block(0, 0, 3, N - 1) = symmetryData().at(shape).coordinates;
  shapeCoordinates.col(N - 1) = Eigen::Vector3d::Zero();
  // Normalize the coordinates
  shapeCoordinates = normalize(shapeCoordinates);

  double permutationalMinimum = std::numeric_limits<double>::max();

  Eigen::Matrix<double, 3, Eigen::Dynamic> permutedShape(3, N);
  std::vector<unsigned> bestPermutation;
  Eigen::Matrix3d bestRotationMatrix;
  do {
    // Construct a permuted shape positions matrix
    for(unsigned i = 0; i < N; ++i) {
      permutedShape.col(permutation.at(i)) = shapeCoordinates.col(i);
    }

    // Perform a quaternion fit
    auto rotationMatrix = fitQuaternion(normalizedPositions, permutedShape);
    permutedShape = rotationMatrix * permutedShape;

    const double value = (normalizedPositions - permutedShape).colwise().squaredNorm().sum();
    if(value < permutationalMinimum) {
      bestPermutation = permutation;
      permutationalMinimum = value;
      bestRotationMatrix = rotationMatrix;
    }
  } while(std::next_permutation(std::begin(permutation), std::end(permutation)));

  for(unsigned i = 0; i < N; ++i) {
    permutedShape.col(bestPermutation.at(i)) = shapeCoordinates.col(i);
  }
  permutedShape = bestRotationMatrix * permutedShape;

  constexpr double scalingLowerBound = 0.5;
  constexpr double scalingUpperBound = 1.1;

  auto scalingMinimizationResult = boost::math::tools::brent_find_minima(
    [&](const double scaling) -> double {
      return (normalizedPositions - scaling * permutedShape).colwise().squaredNorm().sum();
    },
    scalingLowerBound,
    scalingUpperBound,
    std::numeric_limits<double>::digits
  );

  const double normalization = normalizedPositions.colwise().squaredNorm().sum();

  return 100 * scalingMinimizationResult.second / normalization;
}

using NarrowType = std::pair<double, std::unordered_map<unsigned, unsigned>>;

NarrowType shapeHeuristicsNarrow(
  const PositionCollection& stator,
  const PositionCollection& rotor,
  std::unordered_map<unsigned, unsigned> permutation,
  std::vector<unsigned> freeLeftVertices,
  std::vector<unsigned> freeRightVertices
) {
  const unsigned N = stator.cols();

  /* For each left-right mapping, calculate the estimated increased square norm
   * penalty to the rotational fit.
   */
  const unsigned V = freeLeftVertices.size();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> costs (V, V);
  for(unsigned i = 0; i < V; ++i) {
    for(unsigned j = 0; j < V; ++j) {
      costs(i, j) = (
        stator.col(freeLeftVertices.at(i))
        - rotor.col(freeRightVertices.at(j))
      ).squaredNorm();
    }
  }

  /* Generate a mapping that maps free "left" vertices onto the free
   * "right" vertices. Every permutation  of this mapping represents a
   * particular mapping choice whose expected rotational fit penalty we can sum
   * up from the previous matrix.
   */
  auto subPermutation = temple::iota<unsigned>(V);
  decltype(subPermutation) bestPermutation;
  double minimalCost = std::numeric_limits<double>::max();
  do {
    double cost = 0.0;
    for(unsigned i = 0; i < V; ++i) {
      cost += costs(i, subPermutation.at(i));
    }

    if(cost < minimalCost) {
      minimalCost = cost;
      bestPermutation = subPermutation;
    }
  } while(std::next_permutation(std::begin(subPermutation), std::end(subPermutation)));

  // Fuse permutation and best subpermutation
  for(unsigned i = 0; i < V; ++i) {
    permutation.emplace(freeLeftVertices.at(i), freeRightVertices.at(bestPermutation.at(i)));
  }

  assert(permutation.size() == N);

  /* Perform a final fit using all determined pairs and calculate the exact
   * rotational fit penalty
   */
  auto R = fitQuaternion(stator, rotor, permutation);
  auto rotated = R * rotor;

  const double energy = temple::accumulate(
    temple::adaptors::range(N),
    0.0,
    [&](const double carry, const unsigned i) -> double {
      return carry + (
        stator.col(i) - rotated.col(permutation.at(i))
      ).squaredNorm();
    }
  );

  return {energy, permutation};
}

double shapeHeuristics(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
  /* Heuristics employed:
   *
   * - Minimize over the rotation first, then take that solution and minimize
   *   over the scaling factor.
   * - For each tuple of five positions, align the positions, then greedily
   *   choose the best next sequence alignment until all positions are matched.
   *   If there are close seconds in choosing the best next sequence alignment,
   *   branch and explore all of those too.
   *
   * Overall, this works really well for small distortions. It has low deviations
   * from the true error for large deviations, where it suffers from its
   * greediness and not realigning after choosing a sequence alignment.
   *
   * This might also be improved by some pruning criteria (i.e. track the
   * accumulating penalty and the best value found so far and prune if the
   * penalty is significantly higher than the best found).
   *
   * Four positions works reasonably well, but the costs matrix calculated
   * later is not well-converged and the minimal solution is not the shortest
   * cost path through the graph.
   */

  const unsigned N = normalizedPositions.cols();

  if(N < 5) {
    throw std::logic_error("Do not call this heuristics function for less than 5 vertices");
  }

  // Add origin to shape coordinates and renormalize
  PositionCollection shapeCoords (3, N);
  shapeCoords.leftCols(N - 1) = symmetryData().at(shape).coordinates;
  shapeCoords.col(N - 1) = Eigen::Vector3d::Zero();
  shapeCoords = normalize(shapeCoords);

  NarrowType minimalNarrow {std::numeric_limits<double>::max(), {}};

  std::unordered_map<unsigned, unsigned> permutation;
  for(unsigned i = 0; i < N; ++i) {
    permutation[0] = i;
    for(unsigned j = 0; j < N; ++j) {
      if(j == i) {
        continue;
      }

      permutation[1] = j;
      for(unsigned k = 0; k < N; ++k) {
        if(k == i || k == j) {
          continue;
        }

        permutation[2] = k;
        for(unsigned l = 0; l < N; ++l) {
          if(l == i || l == k || l == j) {
            continue;
          }

          permutation[3] = l;

          for(unsigned m = 0; m < N; ++m) {
            if(m == i || m == j || m == k || m == l) {
              continue;
            }

            permutation[4] = m;
            Eigen::Matrix3d R = fitQuaternion(normalizedPositions, shapeCoords, permutation);
            auto rotatedShape = R * shapeCoords;

            /* If the total penalty of a five positions fit is already larger
             * than the tracked minimal penalty, we can discard it already
             * as it can only increase
             */
            double penalty = (
              (normalizedPositions.col(0) - rotatedShape.col(i)).squaredNorm()
              + (normalizedPositions.col(1) - rotatedShape.col(j)).squaredNorm()
              + (normalizedPositions.col(2) - rotatedShape.col(k)).squaredNorm()
              + (normalizedPositions.col(3) - rotatedShape.col(l)).squaredNorm()
              + (normalizedPositions.col(4) - rotatedShape.col(m)).squaredNorm()
            );

            if(penalty > minimalNarrow.first) {
              continue;
            }

            /* Solve the permutational (N-5)! subproblem without realigning all
             * positions.
             */
            std::vector<unsigned> freeLeftVertices;
            freeLeftVertices.reserve(N - 5);
            for(unsigned a = 5; a < N; ++a) {
              freeLeftVertices.push_back(a);
            }
            std::vector<unsigned> freeRightVertices;
            freeRightVertices.reserve(N - 5);
            for(unsigned a = 0; a < N; ++a) {
              if(a != i && a != j && a != k && a != l && a != m) {
                freeRightVertices.push_back(a);
              }
            }

            auto narrowed = shapeHeuristicsNarrow(
              normalizedPositions,
              rotatedShape,
              permutation,
              std::move(freeLeftVertices),
              std::move(freeRightVertices)
            );

            if(narrowed.first < minimalNarrow.first) {
              minimalNarrow = narrowed;
            }
          }
        }
      }
    }
  }

  /* Given the best permutation for the rotational fit, we still have to
   * minimize over the isotropic scaling factor. It is cheaper to reorder the
   * positions once here for the minimization so that memory access is in-order
   * during the repeated scaling minimization function call.
   */

  const auto& bestPermutation = minimalNarrow.second;
  PositionCollection permutedShape(3, N);
  for(unsigned i = 0; i < N; ++i) {
    permutedShape.col(i) = shapeCoords.col(bestPermutation.at(i));
  }
  auto R = fitQuaternion(normalizedPositions, permutedShape);
  permutedShape = R * permutedShape;

  constexpr double scalingLowerBound = 0.5;
  constexpr double scalingUpperBound = 1.1;

  auto scalingMinimizationResult = boost::math::tools::brent_find_minima(
    [&](double scaling) { return (normalizedPositions - scaling * permutedShape).colwise().squaredNorm().sum(); },
    scalingLowerBound,
    scalingUpperBound,
    std::numeric_limits<double>::digits
  );

  const double normalization = normalizedPositions.colwise().squaredNorm().sum();

  return 100 * scalingMinimizationResult.second / normalization;
}

double shape(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
#ifdef NDEBUG
  // In release builds, use heuristics starting from size 9
  constexpr unsigned minSizeForHeuristics = 9;
#else
  // In debug builds, use heuristics starting from size 6
  constexpr unsigned minSizeForHeuristics = 6;
#endif
  if(size(shape) >= minSizeForHeuristics) {
    return shapeHeuristics(normalizedPositions, shape);
  }

  return shapeAlternateImplementation(normalizedPositions, shape);
}

double minimumDistortionAngle(const Shape a, const Shape b) {
  if(size(a) != size(b)) {
    throw std::logic_error("Shapes are not of identical size!");
  }

  // Add origin to shape b's coordinates
  const unsigned S = size(b);
  PositionCollection p (3, S + 1);
  p.block(0, 0, 3, S) = symmetryData().at(b).coordinates;
  p.col(S) = Eigen::Vector3d::Zero();

  return std::asin(
    std::sqrt(shape(normalize(p), a)) / 10
  );
}

double minimalDistortionPathDeviation(
  const PositionCollection& positions,
  Shape a,
  Shape b,
  const double minimumDistortionAngle
) {
  const auto normalized = normalize(positions);
  return (
    std::asin(std::sqrt(shape(normalized, a)) / 10)
    + std::asin(std::sqrt(shape(normalized, b)) / 10)
  ) / minimumDistortionAngle - 1;
}

double minimalDistortionPathDeviation(
  const PositionCollection& positions,
  const Shape a,
  const Shape b
) {
  return minimalDistortionPathDeviation(
    positions,
    a,
    b,
    minimumDistortionAngle(a, b)
  );
}

} // namespace continuous
} // namespace Symmetry
} // namespace Scine
