/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Shapes/ContinuousMeasures.h"

#include "Molassembler/Shapes/Diophantine.h"
#include "Molassembler/Shapes/Partitioner.h"
#include "Molassembler/Shapes/Data.h"

#include "Molassembler/Temple/Adaptors/Iota.h"
#include "Molassembler/Temple/Adaptors/Transform.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Loops.h"
#include "Molassembler/Temple/Optimization/SO3NelderMead.h"
#include "Molassembler/Temple/constexpr/Jsf.h"
#include "Molassembler/Temple/constexpr/Numeric.h"

#include "boost/optional.hpp"
#include "boost/math/tools/minima.hpp"
#include "boost/math/distributions/beta.hpp"
#include <Eigen/Eigenvalues>

#include <random>

namespace Scine {
namespace Molassembler {
namespace Shapes {
namespace Continuous {
namespace {

template<typename PRNG>
Eigen::Vector3d randomVectorOnSphere(const double radius, PRNG& prng) {
  std::normal_distribution<double> normal {};
  Eigen::Vector3d v = Eigen::Vector3d::Zero();

  while(v.norm() < 0.01) {
    v << normal(prng),
         normal(prng),
         normal(prng);
  }

  return radius * v / v.norm();
}

template<typename PRNG>
Eigen::Vector3d randomVectorInSphere(const double radius, PRNG& prng) {
  std::uniform_real_distribution<double> uniform {};
  std::normal_distribution<double> normal {};
  const double u = std::cbrt(uniform(prng));
  Eigen::Vector3d v = Eigen::Vector3d::Zero();

  while(v.norm() < 0.01) {
    v << normal(prng),
         normal(prng),
         normal(prng);
  }

  return radius * u * v / v.norm();
}

template<typename PRNG>
Eigen::Vector3d normallyDistributedVectorInSphere(const double radius, PRNG& prng) {
  std::normal_distribution<double> radiusDistribution {1.0, 0.2};
  std::normal_distribution<double> normal {};
  Eigen::Vector3d v = Eigen::Vector3d::Zero();

  while(v.norm() < 0.01) {
    v << normal(prng),
         normal(prng),
         normal(prng);
  }

  return radius * (1 + radiusDistribution(prng)) * v / v.norm();
}

template<typename PRNG>
Eigen::Matrix<double, 3, Eigen::Dynamic> generateCoordinates(unsigned P, PRNG& prng) {
  Eigen::Matrix<double, 3, Eigen::Dynamic> positions(3, P + 1);
  for(unsigned i = 0; i < P; ++i) {
    positions.col(i) = normallyDistributedVectorInSphere(1.0, prng);
  }
  positions.col(P) = Eigen::Vector3d::Zero();
  return positions;
}

template<typename PRNG>
double sample(const Shapes::Shape shape, PRNG& prng) {
  auto positions = generateCoordinates(Shapes::size(shape), prng);
  positions = Shapes::Continuous::normalize(positions);
  return Shapes::Continuous::shapeCentroidLast(positions, shape).measure;
}

} // namespace


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
  const std::unordered_map<Vertex, Vertex, boost::hash<Vertex>>& p
) {
  assert(centroidIsZero(stator));
  assert(centroidIsZero(rotor));

  Eigen::Matrix4d b = Eigen::Matrix4d::Zero();
  // generate decomposable matrix per atom and add them
  for(const auto& iterPair : p) {
    const auto& statorCol = stator.col(iterPair.first);
    const auto& rotorCol = rotor.col(iterPair.second);

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
    Temple::all_of(
      Temple::iota<unsigned>(P),
      [&](const unsigned i) -> bool {
        return positions.col(i).norm() <= (1 + 1e-5);
      }
    ) && Temple::any_of(
      Temple::iota<unsigned>(P),
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

namespace Fixed {

double element(
  const PositionCollection& normalizedPositions,
  const Elements::Rotation& rotation
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
  if(!Diophantine::has_solution(diophantineConstants, P)) {
    return 100;
  }

  std::vector<unsigned> diophantineMultipliers;
  if(!Diophantine::first_solution(diophantineMultipliers, diophantineConstants, P)) {
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
          auto partitionParticleIndices = Temple::map(partitionIndices, Temple::Functor::at(indicesToPartition));
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
  } while(Diophantine::next_solution(diophantineMultipliers, diophantineConstants, P));

  return 100 * value;
}

double element(
  const PositionCollection& normalizedPositions,
  const Elements::Reflection& reflection
) {
  /* calculateReflectionCSM could be memoized across its a,b arguments, and
   * so could other elements' basic calculation fns. But no reason to optimize
   * since continuous symmetry measures are not a central workload
   */

  const unsigned P = normalizedPositions.cols();
  assert(P > 0);

  Eigen::Hyperplane<double, 3> plane(reflection.normal, 0);
  const Eigen::Matrix3d reflectMatrix = reflection.matrix();
  assert(reflectMatrix.inverse().isApprox(reflectMatrix, 1e-10));

  const std::vector<unsigned> diophantineConstants {2, 1};
  assert(Diophantine::has_solution(diophantineConstants, P));
  std::vector<unsigned> diophantineMultipliers;
  if(!Diophantine::first_solution(diophantineMultipliers, diophantineConstants, P)) {
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
      const double csm = Temple::sum(
        Temple::Adaptors::transform(
          Temple::Adaptors::range(P),
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
  } while(Diophantine::next_solution(diophantineMultipliers, diophantineConstants, P));

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

} // namespace Fixed

std::pair<double, Elements::Rotation> element(
  const PositionCollection& normalizedPositions,
  Elements::Rotation rotation
) {
  using MinimizerType = Temple::SO3NelderMead<>;

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
      return Fixed::element(R * normalizedPositions, rotation);
    },
    NelderMeadChecker {}
  );

  rotation.axis = simplex.at(minimizationResult.minimalIndex).inverse() * rotation.axis;

  return {
    minimizationResult.value,
    rotation
  };
}


std::pair<double, Elements::Reflection> element(
  const PositionCollection& normalizedPositions,
  Elements::Reflection reflection
) {
  using MinimizerType = Temple::SO3NelderMead<>;

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
      return Fixed::element(R * normalizedPositions, reflection);
    },
    NelderMeadChecker {}
  );

  return {
    minimizationResult.value,
    Elements::Reflection {
      simplex.at(minimizationResult.minimalIndex).inverse() * reflection.normal
    }
  };
}

double element(
  const PositionCollection& normalizedPositions,
  const Elements::Inversion& /* inversion */
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
      return Fixed::Cinf(rotatedCoordinates, Eigen::Vector3d::UnitZ());
    }
  };

  Functor functor {normalizedPositions};

  using MinimizerType = Temple::SO3NelderMead<>;

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
  const Elements::ElementGrouping& elementGrouping
) {
  const unsigned p = particles.size();
  const unsigned l = elementGrouping.groups.front().size();
  assert(p * l == unfoldMatrices.cols() / 3);
  assert(Temple::all_of(elementGrouping.groups, [l](const auto& group) { return group.size() == l; }));

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
  const std::vector<Elements::ElementGrouping>& elementGroupings
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
  const Elements::NpGroupingsMapType npGroups;

  OrientationCSMFunctor(
    const PositionCollection& normalizedPositions,
    const PointGroup group
  ) : coordinates(normalizedPositions),
      unfoldMatrices(makeUnfoldMatrices(Elements::symmetryElements(group))),
      foldMatrices(makeFoldMatrices(unfoldMatrices)),
      npGroups(Elements::npGroupings(Elements::symmetryElements(group)))
  {}

  static MatrixType makeUnfoldMatrices(const Elements::ElementsList& elements) {
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
    if(!Diophantine::first_solution(subdivisionMultipliers, subdivisionGroupSizes, P)) {
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
                auto partitionParticles = Temple::map(
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
              auto partitionParticles = Temple::map(
                partitionIndices,
                [&](const unsigned indexOfParticle) -> unsigned {
                  return particleIndices.at(
                    sameSizeParticleIndices.at(indexOfParticle)
                  );
                }
              );

              const double permutationalGroupCSM = groupedSymmetryElements(
                positions,
                std::move(partitionParticles),
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
    } while(Diophantine::next_solution(subdivisionMultipliers, subdivisionGroupSizes, P));

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
    if(Diophantine::has_solution(subdivisionGroupSizes, P)) {
      return diophantine_csm(
        positions,
        subdivisionGroupSizes,
        Temple::iota<unsigned>(P)
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

  const auto elements = Elements::symmetryElements(group);
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

  if(!Diophantine::has_solution(subdivisionGroupSizes, P)) {
    throw std::logic_error("Cannot calculate a CSM for this number of points and this point group. This is most likely an implementation error.");
  }

  OrientationCSMFunctor functor {normalizedPositions, group};

  using MinimizerType = Temple::SO3NelderMead<>;
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

ShapeResult shapeFaithfulPaperImplementation(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
  assert(isNormalized(normalizedPositions));
  const unsigned N = normalizedPositions.cols();

  if(N != size(shape) + 1) {
    throw std::logic_error("Mismatched number of positions between supplied coordinates and shape!");
  }

  auto permutation = Temple::iota<Vertex>(N);

  // Add the origin
  Matrix shapeCoordinates(3, N);
  shapeCoordinates.block(0, 0, 3, N - 1) = coordinates(shape);
  shapeCoordinates.col(N - 1) = Eigen::Vector3d::Zero();
  // Normalize the coordinates
  shapeCoordinates = normalize(shapeCoordinates);

  constexpr double scalingLowerBound = 0.5;
  constexpr double scalingUpperBound = 1.1;

  double permutationalMinimum = std::numeric_limits<double>::max();

  Eigen::Matrix<double, 3, Eigen::Dynamic> permutedShape(3, N);
  std::vector<Vertex> bestPermutation;
  Eigen::Matrix3d bestRotationMatrix;
  do {
    // Construct a permuted shape positions matrix
    for(unsigned i = 0; i < N; ++i) {
      permutedShape.col(i) = shapeCoordinates.col(permutation.at(i));
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

    if(scalingMinimizationResult.second < permutationalMinimum) {
      permutationalMinimum = scalingMinimizationResult.second;
      bestPermutation = permutation;
    }
  } while(std::next_permutation(std::begin(permutation), std::end(permutation)));

  const double normalization = normalizedPositions.colwise().squaredNorm().sum();
  return {
    std::move(bestPermutation),
    100 * permutationalMinimum / normalization
  };
}

namespace Detail {

template<typename EndFunctor>
ShapeResult shapeAlternateImplementationBase(
  const PositionCollection& normalizedPositions,
  const Shape shape,
  EndFunctor&& endFunctor
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

  auto permutation = Temple::iota<Vertex>(N);

  // Add the origin
  Matrix shapeCoordinates(3, N);
  shapeCoordinates.block(0, 0, 3, N - 1) = coordinates(shape);
  shapeCoordinates.col(N - 1) = Eigen::Vector3d::Zero();
  // Normalize the coordinates
  shapeCoordinates = normalize(shapeCoordinates);

  double permutationalMinimum = std::numeric_limits<double>::max();

  Eigen::Matrix<double, 3, Eigen::Dynamic> permutedShape(3, N);
  std::vector<Vertex> bestPermutation;
  Eigen::Matrix3d bestRotationMatrix;
  do {
    // Construct a permuted shape positions matrix
    for(unsigned i = 0; i < N; ++i) {
      permutedShape.col(i) = shapeCoordinates.col(permutation.at(i));
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
  } while(std::next_permutation(std::begin(permutation), endFunctor(permutation)));

  for(unsigned i = 0; i < N; ++i) {
    permutedShape.col(i) = shapeCoordinates.col(bestPermutation.at(i));
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

  return {
    bestPermutation,
    100 * scalingMinimizationResult.second / normalization
  };
}

struct EndIter {
  template<typename T>
  auto operator() (T& t) {
    return std::end(t);
  }
};

struct BackIter {
  template<typename T>
  auto operator() (T& t) {
    assert(!t.empty());
    return --std::end(t);
  }
};

} // namespace Detail

ShapeResult shapeAlternateImplementation(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
  return Detail::shapeAlternateImplementationBase(
    normalizedPositions,
    shape,
    Detail::EndIter {}
  );
}

ShapeResult shapeAlternateImplementationCentroidLast(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
  return Detail::shapeAlternateImplementationBase(
    normalizedPositions,
    shape,
    Detail::BackIter {}
  );
}

using PartialMapping = std::unordered_map<Vertex, Vertex, boost::hash<Vertex>>;
using NarrowType = std::pair<double, PartialMapping>;

NarrowType shapeHeuristicsNarrow(
  const PositionCollection& stator,
  const PositionCollection& rotor,
  std::unordered_map<Vertex, Vertex, boost::hash<Vertex>> permutation,
  std::vector<Vertex> freeLeftVertices,
  std::vector<Vertex> freeRightVertices
) {
  const unsigned N = stator.cols();

  /* For each left-right mapping, calculate the estimated increased square norm
   * penalty to the rotational fit.
   */
  const unsigned V = freeLeftVertices.size();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> costs (V, V);
  for(Vertex i {0}; i < V; ++i) {
    for(Vertex j {0}; j < V; ++j) {
      costs(static_cast<unsigned>(i), static_cast<unsigned>(j)) = (
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
  auto subPermutation = Temple::iota<Vertex>(V);
  decltype(subPermutation) bestPermutation;
  double minimalCost = std::numeric_limits<double>::max();
  do {
    double cost = 0.0;
    for(Vertex i {0}; i < V; ++i) {
      cost += costs(static_cast<unsigned>(i), static_cast<unsigned>(subPermutation.at(i)));
    }

    if(cost < minimalCost) {
      minimalCost = cost;
      bestPermutation = subPermutation;
    }
  } while(std::next_permutation(std::begin(subPermutation), std::end(subPermutation)));

  // Fuse permutation and best subpermutation
  for(Vertex i {0}; i < V; ++i) {
    permutation.emplace(freeLeftVertices.at(i), freeRightVertices.at(bestPermutation.at(i)));
  }

  assert(permutation.size() == N);

  /* Perform a final fit using all determined pairs and calculate the exact
   * rotational fit penalty
   */
  auto R = fitQuaternion(stator, rotor, permutation);
  auto rotated = R * rotor;

  const double energy = Temple::accumulate(
    Temple::Adaptors::range(Vertex(N)),
    0.0,
    [&](const double carry, const Vertex i) -> double {
      return carry + (
        stator.col(i) - rotated.col(permutation.at(i))
      ).squaredNorm();
    }
  );

  return {energy, permutation};
}

ShapeResult shapeHeuristics(
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
   *
   * Further optimization opportunities
   * - Variant with rotation memory (divides number of quaternion fits by the
   *   number of rotations of the iota permutation in the shape)
   */

  const unsigned N = normalizedPositions.cols();

  if(N < 5) {
    throw std::logic_error("Do not call this heuristics function for less than 5 vertices");
  }

  // Add origin to shape coordinates and renormalize
  PositionCollection shapeCoords (3, N);
  shapeCoords.leftCols(N - 1) = coordinates(shape);
  shapeCoords.col(N - 1) = Eigen::Vector3d::Zero();
  shapeCoords = normalize(shapeCoords);

  NarrowType minimalNarrow {std::numeric_limits<double>::max(), {}};
  PartialMapping permutation;

  // i != j != k != l != m with {i, j, k, l, m} in [0, N)
  Temple::Loops::different(
    [&](const std::vector<Vertex>& vertices) {
      permutation.clear();
      for(unsigned i = 0; i < 5; ++i) {
        permutation.emplace(i, vertices[i]);
      }

      const Vertex i = vertices[0];
      const Vertex j = vertices[1];
      const Vertex k = vertices[2];
      const Vertex l = vertices[3];
      const Vertex m = vertices[4];

      Eigen::Matrix3d R = fitQuaternion(normalizedPositions, shapeCoords, permutation);
      auto rotatedShape = R * shapeCoords;

      /* If the total penalty of a five positions fit is already larger
       * than the tracked minimal penalty, we can discard it already
       * as it can only increase
       */
      const double penalty = (
        (normalizedPositions.col(0) - rotatedShape.col(i)).squaredNorm()
        + (normalizedPositions.col(1) - rotatedShape.col(j)).squaredNorm()
        + (normalizedPositions.col(2) - rotatedShape.col(k)).squaredNorm()
        + (normalizedPositions.col(3) - rotatedShape.col(l)).squaredNorm()
        + (normalizedPositions.col(4) - rotatedShape.col(m)).squaredNorm()
      );

      if(penalty > minimalNarrow.first) {
        return;
      }

      /* Solve the permutational (N-5)! subproblem without realigning all
       * positions.
       */
      std::vector<Vertex> freeLeftVertices;
      freeLeftVertices.reserve(N - 5);
      for(Vertex a {5}; a < N; ++a) {
        freeLeftVertices.push_back(a);
      }
      std::vector<Vertex> freeRightVertices;
      freeRightVertices.reserve(N - 5);
      for(Vertex a {0}; a < N; ++a) {
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
    },
    5,
    Vertex {N}
  );

  /* Given the best permutation for the rotational fit, we still have to
   * minimize over the isotropic scaling factor. It is cheaper to reorder the
   * positions once here for the minimization so that memory access is in-order
   * during the repeated scaling minimization function call.
   */

  std::vector<Vertex> bestPermutation(minimalNarrow.second.size());
  for(const auto& iterPair : minimalNarrow.second) {
    bestPermutation.at(iterPair.first) = iterPair.second;
  }

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

  return {
    std::move(bestPermutation),
    100 * scalingMinimizationResult.second / normalization
  };
}

ShapeResult shapeHeuristicsCentroidLast(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
  /* Same as shapeHeuristics, except we know two things:
   *
   * The centroid mapping is known, so the centroid index is unavailable in
   * mapping. Additionally, we have a definitive mapping that we can exploit in
   * every five-index quaternion fit. So the complexity reduces to
   * (N - 1)! / (N - 5)! quaternion fits, I think.
   */

  const unsigned N = normalizedPositions.cols();

  if(N < 6) {
    throw std::logic_error("Do not call this heuristics function for less than 5 vertices");
  }

  // Add origin to shape coordinates and renormalize
  PositionCollection shapeCoords (3, N);
  shapeCoords.leftCols(N - 1) = coordinates(shape);
  shapeCoords.col(N - 1) = Eigen::Vector3d::Zero();
  shapeCoords = normalize(shapeCoords);

  NarrowType minimalNarrow {std::numeric_limits<double>::max(), {}};
  PartialMapping permutation;

  // i != j != k != l != m with {i, j, k, l, m} in [0, N - 1)
  Temple::Loops::different(
    [&](const std::vector<Vertex>& vertices) {
      permutation.clear();
      permutation.emplace(N - 1, N - 1);
      for(unsigned i = 0; i < 5; ++i) {
        permutation.emplace(i, vertices[i]);
      }

      const Vertex i = vertices[0];
      const Vertex j = vertices[1];
      const Vertex k = vertices[2];
      const Vertex l = vertices[3];
      const Vertex m = vertices[4];

      Eigen::Matrix3d R = fitQuaternion(normalizedPositions, shapeCoords, permutation);
      auto rotatedShape = R * shapeCoords;

      /* If the total penalty of a five positions fit is already larger
       * than the tracked minimal penalty, we can discard it already
       * as it can only increase
       */
      const double penalty = (
        (normalizedPositions.col(0) - rotatedShape.col(i)).squaredNorm()
        + (normalizedPositions.col(1) - rotatedShape.col(j)).squaredNorm()
        + (normalizedPositions.col(2) - rotatedShape.col(k)).squaredNorm()
        + (normalizedPositions.col(3) - rotatedShape.col(l)).squaredNorm()
        + (normalizedPositions.col(4) - rotatedShape.col(m)).squaredNorm()
        + (normalizedPositions.col(N - 1) - rotatedShape.col(N - 1)).squaredNorm()
      );

      if(penalty > minimalNarrow.first) {
        return;
      }

      /* Solve the permutational (N-5)! subproblem without realigning all
       * positions.
       */
      std::vector<Vertex> freeLeftVertices;
      freeLeftVertices.reserve(N - 6);
      for(Vertex a {5}; a < N - 1; ++a) {
        freeLeftVertices.push_back(a);
      }
      std::vector<Vertex> freeRightVertices;
      freeRightVertices.reserve(N - 6);
      for(Vertex a {0}; a < N - 1; ++a) {
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
    },
    5,
    Vertex {N - 1}
  );

  /* Given the best permutation for the rotational fit, we still have to
   * minimize over the isotropic scaling factor. It is cheaper to reorder the
   * positions once here for the minimization so that memory access is in-order
   * during the repeated scaling minimization function call.
   */

  std::vector<Vertex> bestPermutation(minimalNarrow.second.size());
  for(const auto& iterPair : minimalNarrow.second) {
    bestPermutation.at(iterPair.first) = iterPair.second;
  }

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

  return {
    std::move(bestPermutation),
    100 * scalingMinimizationResult.second / normalization
  };
}


ShapeResult shape(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
#ifdef NDEBUG
  // In release builds, use heuristics starting from size 9
  constexpr unsigned minSizeForHeuristics = 8;
#else
  // In debug builds, use heuristics starting from size 6
  constexpr unsigned minSizeForHeuristics = 6;
#endif
  if(size(shape) >= minSizeForHeuristics) {
    return shapeHeuristics(normalizedPositions, shape);
  }

  return shapeAlternateImplementation(normalizedPositions, shape);
}

ShapeResult shapeCentroidLast(
  const PositionCollection& normalizedPositions,
  const Shape shape
) {
#ifdef NDEBUG
  // In release builds, use heuristics starting from size 8
  constexpr unsigned minSizeForHeuristics = 8;
#else
  // In debug builds, use heuristics starting from size 6
  constexpr unsigned minSizeForHeuristics = 6;
#endif
  if(size(shape) >= minSizeForHeuristics) {
    return shapeHeuristicsCentroidLast(normalizedPositions, shape);
  }

  return shapeAlternateImplementationCentroidLast(normalizedPositions, shape);
}

double minimumDistortionAngle(const Shape a, const Shape b) {
  if(size(a) != size(b)) {
    throw std::logic_error("Shapes are not of identical size!");
  }

  // Add origin to shape b's coordinates
  const unsigned S = size(b);
  PositionCollection p (3, S + 1);
  p.block(0, 0, 3, S) = coordinates(b);
  p.col(S) = Eigen::Vector3d::Zero();

  return std::asin(
    std::sqrt(shape(normalize(p), a).measure) / 10
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
    std::asin(std::sqrt(shape(normalized, a).measure) / 10)
    + std::asin(std::sqrt(shape(normalized, b).measure) / 10)
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

std::array<double, 4> randomCloudDistributionParameters(
  const Shape shape,
  const unsigned N,
  const unsigned seed
) {
  Temple::JSF64 prng;
  prng.seed(seed);

  // Generate samples
  std::vector<double> samples(N);

  for(unsigned i = 0; i < N; ++i) {
    samples[i] = sample(shape, prng);
  }

  // Fit a distribution
  const double scale = Temple::max(samples);
  const double samplesMean = Temple::average(samples);
  const double samplesStddev = Temple::stddev(samples, samplesMean);

  // Scale to [0, 1] for beta distribution work
  const double betaMean = samplesMean / scale;
  const double betaStddev = samplesStddev / scale;
  const double betaVariance = betaStddev * betaStddev;

  using Beta = boost::math::beta_distribution<>;
  const double a = Beta::find_alpha(betaMean, betaVariance);
  const double b = Beta::find_beta(betaMean, betaVariance);

  return {{a, b, 0.0, scale}};
}

boost::optional<double> probabilityRandomCloud(const double measure, const Shape shape) {
  // Generated with N = 200
  static const std::unordered_map<Shapes::Shape, std::array<double, 4>> randomDistributionParameters {
    {Shapes::Shape::Line, {{0.766124, 1.91041, 0, 125.844}}},
    {Shapes::Shape::Bent, {{0.29139, 1.89203, 0, 80.4878}}},
    {Shapes::Shape::EquilateralTriangle, {{1.29541, 3.38419, 0, 96.0021}}},
    {Shapes::Shape::VacantTetrahedron, {{1.24833, 3.86221, 0, 87.781}}},
    {Shapes::Shape::T, {{1.04193, 3.82784, 0, 78.2735}}},
    {Shapes::Shape::Tetrahedron, {{3.72178, 6.7318, 0, 92.1706}}},
    {Shapes::Shape::Square, {{2.86648, 4.12431, 0, 65.1664}}},
    {Shapes::Shape::Seesaw, {{2.46897, 6.3599, 0, 73.4791}}},
    {Shapes::Shape::TrigonalPyramid, {{2.2215, 3.2188, 0, 60.5622}}},
    {Shapes::Shape::SquarePyramid, {{4.328, 5.55626, 0, 55.414}}},
    {Shapes::Shape::TrigonalBipyramid, {{3.68118, 4.21196, 0, 55.8912}}},
    {Shapes::Shape::Pentagon, {{3.89701, 4.22035, 0, 58.1819}}},
    {Shapes::Shape::Octahedron, {{5.70985, 6.08247, 0, 63.8315}}},
    {Shapes::Shape::TrigonalPrism, {{4.39501, 4.36312, 0, 50.0585}}},
    {Shapes::Shape::PentagonalPyramid, {{5.073, 4.38382, 0, 44.3454}}},
    {Shapes::Shape::Hexagon, {{5.97833, 8.70608, 0, 70.426}}},
    {Shapes::Shape::PentagonalBipyramid, {{5.62079, 4.51257, 0, 47.6727}}},
    {Shapes::Shape::CappedOctahedron, {{5.3341, 6.54619, 0, 55.0139}}},
    {Shapes::Shape::CappedTrigonalPrism, {{3.62137, 3.57516, 0, 48.2583}}},
    {Shapes::Shape::SquareAntiprism, {{6.49925, 5.6291, 0, 45.6774}}},
    {Shapes::Shape::Cube, {{7.38926, 6.12884, 0, 48.6829}}},
    {Shapes::Shape::TrigonalDodecahedron, {{6.44829, 5.57766, 0, 45.2966}}},
    {Shapes::Shape::HexagonalBipyramid, {{6.4292, 4.48488, 0, 42.8605}}},
  };

  const auto findIter = randomDistributionParameters.find(shape);
  if(findIter == std::end(randomDistributionParameters)) {
    return boost::none;
  }

  const double a = findIter->second.at(0);
  const double b = findIter->second.at(1);
  const double scale = findIter->second.at(3);
  if(measure > scale) {
    return 1.0;
  }

  return boost::math::cdf(boost::math::beta_distribution<>(a, b), measure / scale);
}

} // namespace Continuous
} // namespace Shapes
} // namespace Molassembler
} // namespace Scine
