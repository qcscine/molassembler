/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "boost/program_options.hpp"

#include "molassembler/Shapes/ContinuousMeasures.h"

#include "molassembler/Shapes/Data.h"
#include "molassembler/Shapes/PropertyCaching.h"

#include "boost/math/tools/minima.hpp"
#include "molassembler/Temple/constexpr/Numeric.h"
#include "molassembler/Temple/Adaptors/AllPairs.h"
#include "molassembler/Temple/Adaptors/Iota.h"
#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/Stringify.h"
#include "molassembler/Temple/Random.h"
#include "molassembler/Temple/constexpr/Jsf.h"
#include "molassembler/Temple/Permutations.h"
#include "molassembler/Temple/Loops.h"

#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <chrono>

/* Summary
 *
 * This file contains a number of methods trying to cheapen finding the correct
 * vertex mapping to evaluate the continuous shape measure.
 *
 * Note that the problem does have some rotational redundancy depending on the
 * number of rotations that the shape has. But exploiting that gets you maybe
 * (N - 1)!. Crasser solutions are necessary.
 */

using namespace Scine;
using namespace Shapes;
using namespace Continuous;

constexpr unsigned factorial(unsigned x) {
  if(x <= 1) {
    return 1;
  }

  return x * factorial(x - 1);
}

template<std::size_t ... Inds>
constexpr auto makeFactorials(std::index_sequence<Inds ...> /* inds */) {
  return std::array<unsigned, sizeof...(Inds)> {
    factorial(Inds + 1)...
  };
}

constexpr auto facs = makeFactorials(std::make_index_sequence<14> {});

Continuous::PositionCollection addOrigin(const Continuous::PositionCollection& vs) {
  const unsigned N = vs.cols();
  Continuous::PositionCollection positions(3, N + 1);
  for(unsigned i = 0; i < N; ++i) {
    positions.col(i) = vs.col(i);
  }

  // Add origin point explicitly to consideration
  positions.col(N) = Eigen::Vector3d::Zero(3);
  return positions;
}

void distort(Eigen::Ref<Continuous::PositionCollection> positions, const double distortionNorm = 0.01) {
  const unsigned N = positions.cols();
  for(unsigned i = 0; i < N; ++i) {
    positions.col(i) += distortionNorm * Eigen::Vector3d::Random().normalized();
  }
}

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

using Parameters = std::vector<Vertex>;

struct Energy {
  const PositionCollection& referencePositions;
  const PositionCollection& shapePositions;

  double operator() (const Parameters& p) const {
    const unsigned N = shapePositions.cols();

    static PositionCollection permuted(3, N);
    for(unsigned i = 0; i < N; ++i) {
      permuted.col(i) = shapePositions.col(p[i]);
    }

    auto rotationMatrix = fitQuaternion(referencePositions, permuted);
    permuted = rotationMatrix * permuted;

    return (referencePositions - permuted).colwise().squaredNorm().sum();

    // // Minimize over isotropic scaling factor
    // auto scalingMinimizationResult = boost::math::tools::brent_find_minima(
    //   [&](const double scaling) -> double {
    //     return (referencePositions - scaling * permuted).colwise().squaredNorm().sum();
    //   },
    //   0.5,
    //   1.1,
    //   std::numeric_limits<double>::digits
    // );

    // return scalingMinimizationResult.second;
  }
};

struct OldShapeAlgorithm {
  virtual ~OldShapeAlgorithm() = default;

  static PositionCollection shapeCoordinates(const Shape shape) {
    return normalize(addOrigin(coordinates(shape)));
  }

  virtual double shape(const PositionCollection& positions, Shape shape) = 0;
  virtual std::string name() const = 0;
};

/**
 * @brief Abstract base class for algorithms to conform to
 */
struct ShapeAlgorithm {
  virtual ~ShapeAlgorithm() = default;

  using Permutation = std::vector<Vertex>;
  using ResultPair = std::pair<double, Permutation>;

  static Permutation invert(const Permutation& p) {
    const unsigned P = p.size();
    Permutation inverse(P);
    for(unsigned i = 0; i < P; ++i) {
      inverse.at(p.at(i)) = i;
    }
    return inverse;
  }

  /* A forwards permutation maps from positions to shape indices,
   * a backwards permutation maps from shape indices to positions
   *
   * Only backwards permutations have the (added) centroid at the back, and
   * hence their [O, S] = [0, N - 1] = [0, N) ranges can be rotated.
   */
  static std::set<ShapeAlgorithm::Permutation> generateRotations(
    const Permutation& forwardsPermutation,
    const Shape shape
  ) {
    auto backwards = ShapeAlgorithm::invert(forwardsPermutation);
    Vertex centroid = backwards.back();
    backwards.pop_back();
    auto rotations = Properties::generateAllRotations(shape, backwards);

    std::set<ShapeAlgorithm::Permutation> forwardsPermutationRotations;
    for(auto rotation : rotations) {
      rotation.push_back(centroid);
      forwardsPermutationRotations.emplace(
        ShapeAlgorithm::invert(rotation)
      );
    }
    return forwardsPermutationRotations;
  }

  static PositionCollection shapeCoordinates(const Shape shape) {
    return normalize(addOrigin(coordinates(shape)));
  }

  virtual ResultPair shape(
    const PositionCollection& positions,
    Shape shape,
    const Permutation& correctPermutation
  ) = 0;
  virtual std::string name() const = 0;
};

/**
 * @brief Reference continuous shape measure calculation method
 *
 * Complexity: Theta(N!)
 */
struct AllPermutations {
  ShapeAlgorithm::ResultPair shape(const PositionCollection& positions, Shape shape) {
    Energy energy {positions, ShapeAlgorithm::shapeCoordinates(shape)};
    ShapeAlgorithm::Permutation permutation = Temple::iota<Vertex>(Vertex(positions.cols()));
    ShapeAlgorithm::ResultPair minimal {std::numeric_limits<double>::max(), {}};
    do {
      double permutationEnergy = energy(permutation);

      if(permutationEnergy < minimal.first) {
        minimal = {permutationEnergy, permutation};
      }
    } while(Temple::InPlace::next_permutation(permutation));

    return minimal;
  }

  std::string name() const {
    return "Reference";
  }
};

struct Cooling {
  // This works best for annealing and tunneling
  static double linear(const unsigned step, const unsigned steps) {
    return 1.0 - (static_cast<double>(step) / steps);
  }

  // This doesn't work well at all for annealing or tunneling, but it's here
  static double exponential(const unsigned step, const unsigned steps) {
    const double alpha = std::exp(std::log(1e-2) / steps);
    return std::pow(alpha, step);
  }
};

/**
 * @brief Simulated annealing method
 *
 * With well tuned temperature, this works pretty well, but it's stochastic in
 * nature and its good results are due to tracking the lowest energy state
 * visited instead of using the final annealing state. The energy surface for
 * small distortion norms is super jagged and the "correct" state has no
 * low-energy neighbors, making this type of approach problematic. Behaves
 * better for high distortion norms as there state space can be better sampled.
 *
 * Complexity: Needs at least O((N - 2)!) steps to explore enough state space
 * to encounter the minimum.
 */
struct Anneal final : OldShapeAlgorithm {
  using Permutation = std::vector<Vertex>;

  Temple::JSF64 prng;
  std::ofstream trace {"anneal_trace.csv"};

  Anneal() {
    prng.seed(0);
  }

  static double acceptanceProbability(
    const double energy,
    const double prospectiveEnergy,
    const double temperature
  ) {
    if(prospectiveEnergy < energy) {
      return 1.0;
    }

    return std::exp(-(prospectiveEnergy - energy) / temperature);
  }

  static double temperature(const unsigned step, const unsigned steps) {
    return Cooling::linear(step, steps);
  }

  //! Assumes temperature in 0, 1
  Permutation generateMove(Permutation state, const double temperature) {
    const unsigned N = state.size();

    const unsigned nSwaps = 1 + std::round(temperature * Temple::Random::getSingle<unsigned>(0, N - 1, prng));
    for(unsigned s = 0; s < nSwaps; ++s) {
      const unsigned i = Temple::Random::getSingle<unsigned>(0, N - 2, prng);
      std::swap(state.at(i), state.at(i + 1));
    }

    return state;
  }

  double shape(const PositionCollection& positions, Shape shape) final {
    const unsigned steps = facs.at(positions.cols() - 2);
    Energy energy {positions, shapeCoordinates(shape)};

    auto parameters = Temple::iota<Vertex>(positions.cols());
    decltype(parameters) prospectiveParameters;

    double temperatureMultiplier = Temple::accumulate(
      Temple::Adaptors::range(10),
      0.0,
      [&](const double carry, unsigned /* i */) -> double {
        Temple::Random::shuffle(parameters, prng);
        return carry + energy(parameters);
      }
    ) / 10;

    double currentEnergy = energy(parameters);
    double minimalEnergy = currentEnergy;
    for(unsigned step = 0; step < steps; ++step) {
      const double currentTemperature = temperature(step, steps);
      prospectiveParameters = generateMove(parameters, currentTemperature);
      const double prospectiveEnergy = energy(prospectiveParameters);
      const double p = acceptanceProbability(currentEnergy, prospectiveEnergy, temperatureMultiplier * currentTemperature);

      trace << step << ", " << currentEnergy << ", " << temperatureMultiplier * currentTemperature << ", " << prospectiveEnergy << ", " << p << "\n";

      if(p >= Temple::Random::getSingle(0.0, 1.0, prng)) {
        std::swap(parameters, prospectiveParameters);
        currentEnergy = prospectiveEnergy;
      }

      minimalEnergy = std::min(minimalEnergy, prospectiveEnergy);
    }

    return minimalEnergy;
  }

  std::string name() const final {
    return "Annealing";
  }
};

/**
 * @brief Stochastic tunneling method
 *
 * This works by transforming the energy function to level out minima that have
 * higher function values than the lowest found. This then should help losing
 * time in unhelpful minima.
 *
 * Works pretty well even when the temperature is not well tuned. You have to
 * adjust the gamma parameter for the particular problem you have, and the
 * value we have now works well for CShMs.
 *
 * But now that simulated annealing is well tuned and since there aren't really
 * and multi-state minima to get stuck in in this energy surface for small
 * distortions, this doesn't perform any better.
 *
 * The move generator for this and for simulated annealing generate multi-swap
 * moves for high temperatures and increase locality as the temperature is
 * lowered, generating only single adjacent swap moves.
 *
 * Complexity: Needs at least O((N - 2)!) steps to explore enough state space
 * to encounter the minimum.
 */
struct Tunnel final : OldShapeAlgorithm {
  using Permutation = std::vector<Vertex>;

  Temple::JSF64 prng;
  std::ofstream trace {"tunnel_trace.csv"};
  double gamma = 0.1;

  Tunnel() {
    prng.seed(0);
  }

  static double acceptanceProbability(
    const double energy,
    const double prospectiveEnergy,
    const double temperature
  ) {
    if(prospectiveEnergy < energy) {
      return 1.0;
    }

    return std::exp(-(prospectiveEnergy - energy) / temperature);
  }

  //! Assumes parameters are untransformed!
  double energyTransformation(const double energy, const double lowestEnergy) const {
    return 1.0 - std::exp(-gamma * (energy - lowestEnergy));
  }

  static double temperature(const unsigned step, const unsigned steps) {
    return Cooling::linear(step, steps);
  }

  //! Assumes temperature in 0, 1
  Permutation generateMove(Permutation state, const double temperature) {
    const unsigned N = state.size();

    const unsigned nSwaps = 1 + std::round(temperature * Temple::Random::getSingle<unsigned>(0, N - 1, prng));
    for(unsigned s = 0; s < nSwaps; ++s) {
      const unsigned i = Temple::Random::getSingle<unsigned>(0, N - 2, prng);
      std::swap(state.at(i), state.at(i + 1));
    }

    return state;
  }

  double shape(const PositionCollection& positions, Shape shape) final {
    const unsigned steps = facs.at(positions.cols() - 2);
    Energy energy {positions, shapeCoordinates(shape)};

    auto parameters = Temple::iota<Vertex>(positions.cols());
    decltype(parameters) prospectiveParameters;

    double currentEnergy = energy(parameters);
    double lowestEnergy = currentEnergy;
    double transformedCurrentEnergy = energyTransformation(currentEnergy, lowestEnergy);

    for(unsigned step = 0; step < steps; ++step) {
      const double currentTemperature = temperature(step, steps);
      prospectiveParameters = generateMove(parameters, currentTemperature);
      const double prospectiveEnergy = energy(prospectiveParameters);
      const double transformedProspectiveEnergy = energyTransformation(prospectiveEnergy, lowestEnergy);
      const double p = acceptanceProbability(transformedCurrentEnergy, transformedProspectiveEnergy, 0.7 * currentTemperature);

      trace << step << ", " << currentEnergy << ", " << 0.7 * currentTemperature << ", " << prospectiveEnergy << ", " << p << "\n";

      if(p >= Temple::Random::getSingle(0.0, 1.0, prng)) {
        std::swap(parameters, prospectiveParameters);
        currentEnergy = prospectiveEnergy;
        transformedCurrentEnergy = transformedProspectiveEnergy;
      }

      if(prospectiveEnergy < lowestEnergy) {
        lowestEnergy = prospectiveEnergy;
        transformedCurrentEnergy = energyTransformation(currentEnergy, lowestEnergy);
      }
    }

    return lowestEnergy;
  }

  std::string name() const final {
    return "Tunneling";
  }
};

template<typename T, std::size_t N>
struct CircularBuffer {
  void insert(T value) {
    if(size_ < N) {
      buffer_[size_] = std::move(value);
      ++size_;
    } else {
      buffer_[start_] = std::move(value);
      start_ = (start_ + 1) % N;
    }
  }

  T min() const {
    assert(size_ > 0);
    return *std::min_element(
      std::begin(buffer_),
      std::end(buffer_)
    );
  }

  double average() const {
    return static_cast<double>(
      std::accumulate(
        std::begin(buffer_),
        std::end(buffer_),
        T {0},
        std::plus<>()
      )
    ) / size_;
  }

  double variance(const double average) const {
    return static_cast<double>(
      std::accumulate(
        std::begin(buffer_),
        std::end(buffer_),
        T {0},
        [&](const double carry, const double value) -> double {
          return carry + std::pow(value - average, 2);
        }
      )
    ) / size_;
  }

  void clear() {
    size_ = 0;
    start_ = 0;
  }

  std::size_t size() const {
    return size_;
  }

  std::array<T, N> buffer_;
  std::size_t start_ = 0;
  std::size_t size_ = 0;
};

/**
 * @brief An attempt at thermodynamic simulated annealing
 *
 * In thermodynamic simulated annealing, you try to exploit statistical
 * mechanics / thermodynamics to drive your cooling schedule optimally.
 *
 * However, I couldn't ever get this to work. The temperature adjustments are
 * always too large and no amount of messing about ever got them right.
 *
 * Complexity: ???
 */
struct ThermodynamicAnneal final : OldShapeAlgorithm {
  Temple::JSF64 prng;
  using Permutation = std::vector<Vertex>;

  std::unordered_map<unsigned, unsigned> stateIndexReduction;
  std::vector<double> energies;
  Eigen::SparseMatrix<unsigned> Q;
  CircularBuffer<double, 1000> lastEnergies;
  std::ofstream trace {"thermo_trace.csv"};

  ThermodynamicAnneal() {
    prng.seed(0);
  }

  static double acceptanceProbability(
    const double energy,
    const double prospectiveEnergy,
    const double temperature
  ) {
    if(prospectiveEnergy < energy) {
      return 1.0;
    }

    return std::exp(-(prospectiveEnergy - energy) / temperature);
  }

  //! Assumes temperature in 0, 100
  Permutation generateMove(Permutation state, const double temperature) {
    const unsigned N = state.size();

    const unsigned nSwaps = 1 + std::round((temperature / 10) * Temple::Random::getSingle<unsigned>(0, N - 1, prng));
    for(unsigned s = 0; s < nSwaps; ++s) {
      const unsigned i = Temple::Random::getSingle<unsigned>(0, N - 2, prng);
      std::swap(state.at(i), state.at(i + 1));
    }

    return state;
  }

  double partitionFunction(const double temperature) const {
    return Temple::accumulate(
      energies,
      0.0,
      [&](const double carry, const double energy) -> double {
        return carry + std::exp(-energy / temperature);
      }
    );
  }

  template<typename F>
  static double centralDifference(F&& f, const double x, const double h) {
    return (
      f(x + h) - f(x - h)
    ) / (2 * h);
  }

  double averageEnergy(const double temperature) const {
    // Estimate dln Z/dT by central finite difference
    return temperature * temperature * centralDifference(
      [&](const double x) { return std::log(partitionFunction(x)); },
      temperature,
      1e-4
    );
  }

  double heatCapacity(const double temperature) const {
    return centralDifference(
      [&](const double x) { return averageEnergy(x); },
      temperature,
      1e-4
    );
  }

  double relaxationTime(const double temperature) const {
    // Q is strictly lower triangular (i.e. always access i > j)
    Eigen::SparseMatrix<double> G(Q.rows(), Q.cols());
    for(Eigen::Index k = 0; k < Q.outerSize(); ++k) {
      // Sum up all values in Q's column
      double sum = 0;
      for(Eigen::SparseMatrix<unsigned>::InnerIterator it(Q,k); it; ++it) {
        sum += it.value();
      }

      for(Eigen::SparseMatrix<unsigned>::InnerIterator it(Q,k); it; ++it) {
        if(it.row() != it.col()) {
          const double energyDiff = energies.at(it.row()) - energies.at(it.col());
          if(energyDiff > 0) {
            G.coeffRef(it.row(), it.col()) = it.value() * std::exp(-energyDiff / temperature);
          } else {
            G.coeffRef(it.row(), it.col()) = it.value();
          }
        }
      }

      // Diagonal entries of G are 1 - other column's entries
      sum = 0;
      for(Eigen::SparseMatrix<double>::InnerIterator it(G,k); it; ++it) {
        sum += it.value();
      }
      G.coeffRef(k, k) = 1 - sum;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<double>> solver(G);
    return - 1.0 / solver.eigenvalues()(1);
  }

  double updateTemperature(const double temperature, double minimalEnergy) {
    const double hundredEnergiesAverage = lastEnergies.average();
    const double hundredEnergiesVariance = lastEnergies.variance(hundredEnergiesAverage);
    const double thermodynamicSpeed = (hundredEnergiesAverage - minimalEnergy) / hundredEnergiesVariance;

    const double epsilon = relaxationTime(temperature);
    const double C = heatCapacity(temperature);

    double Theta = 1 + temperature * centralDifference([&](double x) { return heatCapacity(x);}, temperature, 1e-4) / (2 * C);
    auto calculateDelta = [&](const double delta) {
      double correction = std::sqrt(1 + (Theta * epsilon * delta) / temperature);
      return - thermodynamicSpeed * temperature / (
        epsilon * std::sqrt(C) * correction
      );
    };

    // NOTE: this can diverge too :(
    double delta = -1.0;
    for(unsigned i = 0; i < 10; ++i) {
      delta = calculateDelta(delta);
    }

    // NOTE 1e-5 is abitrary, the updates are always too large!
    return std::max(0.0, temperature + 1e-5 * delta);
  }

  double shape(const PositionCollection& positions, Shape shape) final {
    const unsigned steps = 5e5;

    stateIndexReduction.clear();
    energies.clear();
    energies.reserve(100);
    Q.resize(0,0);
    Q.reserve(100);
    lastEnergies.clear();

    const unsigned N = positions.cols();
    Energy energy {positions, shapeCoordinates(shape)};

    auto permutation = Temple::iota<Vertex>(N);
    decltype(permutation) prospectiveMove;
    double currentEnergy = energy(permutation);
    double minimalEnergy = currentEnergy;

    unsigned iop = Temple::permutationIndex(permutation);
    stateIndexReduction.emplace(iop, 0);
    energies.push_back(currentEnergy);
    unsigned currentStateIndex = 0;

    double temperature = 10.0;

    /* First 1000 steps of annealing without changing temperature to collect
     * statistics we can use to update the temperature
     */
    unsigned step = 0;
    for(; step < 1000; ++step) {
      prospectiveMove = generateMove(permutation, temperature);
      const double prospectiveEnergy = energy(prospectiveMove);
      const double p = acceptanceProbability(currentEnergy, prospectiveEnergy, temperature);

      iop = Temple::permutationIndex(prospectiveMove);
      unsigned prospectiveStateIndex;
      auto stateIndexFindIter = stateIndexReduction.find(iop);
      if(stateIndexFindIter == std::end(stateIndexReduction)) {
        prospectiveStateIndex = stateIndexReduction.size();
        stateIndexReduction.emplace(iop, prospectiveStateIndex);
        energies.push_back(prospectiveEnergy);
        Q.conservativeResize(energies.size(), energies.size());
        Q.insert(prospectiveStateIndex, currentStateIndex) = 1;
      } else {
        prospectiveStateIndex = stateIndexFindIter->second;
        unsigned col = prospectiveStateIndex;
        unsigned row = currentStateIndex;
        if(col > row) {
          std::swap(col, row);
        }
        Q.coeffRef(row, col) = Q.coeff(row, col) + 1;
      }

      trace << step << ", " << currentEnergy << ", " << temperature << ", " << prospectiveEnergy << ", " << p << "\n";

      if(p >= Temple::Random::getSingle(0.0, 1.0, prng)) {
        std::swap(permutation, prospectiveMove);
        currentEnergy = prospectiveEnergy;
        currentStateIndex = prospectiveStateIndex;
      }

      minimalEnergy = std::min(minimalEnergy, prospectiveEnergy);
      lastEnergies.insert(currentEnergy);
    }

    /* Now let temperature freely vary according to update formula */
    while(temperature > 0.01 && step < steps) {
      temperature = updateTemperature(temperature, minimalEnergy);
      prospectiveMove = generateMove(permutation, temperature);
      const double prospectiveEnergy = energy(prospectiveMove);
      const double p = acceptanceProbability(currentEnergy, prospectiveEnergy, temperature);

      /* Update tracking state */
      iop = Temple::permutationIndex(prospectiveMove);
      unsigned prospectiveStateIndex;
      auto stateIndexFindIter = stateIndexReduction.find(iop);
      if(stateIndexFindIter == std::end(stateIndexReduction)) {
        prospectiveStateIndex = stateIndexReduction.size();
        stateIndexReduction.emplace(iop, prospectiveStateIndex);
        energies.push_back(prospectiveEnergy);
        Q.conservativeResize(energies.size(), energies.size());
        Q.insert(prospectiveStateIndex, currentStateIndex) = 1;
      } else {
        prospectiveStateIndex = stateIndexFindIter->second;
        unsigned col = prospectiveStateIndex;
        unsigned row = currentStateIndex;
        if(col > row) {
          std::swap(col, row);
        }
        Q.coeffRef(row, col) = Q.coeff(row, col) + 1;
      }

      trace << step << ", " << currentEnergy << ", " << temperature << ", " << prospectiveEnergy << ", " << p << "\n";

      // Conditionally accept the move
      if(p >= Temple::Random::getSingle(0.0, 1.0, prng)) {
        std::swap(permutation, prospectiveMove);
        currentEnergy = prospectiveEnergy;
        currentStateIndex = prospectiveStateIndex;
      }

      minimalEnergy = std::min(minimalEnergy, prospectiveEnergy);
      lastEnergies.insert(currentEnergy);
      ++step;
    }

    return minimalEnergy;
  }

  std::string name() const final {
    return "Thermodyn.";
  }
};

/**
 * @brief Fixed-step number greedy minimization with shuffling
 *
 * This just greedily minimizes the permutation and shuffles when it's at the
 * minimum. Unreliable.
 *
 * Complexity: Needs at least Theta(N - 2)! steps to discover enough state
 * space that it might find the minimum.
 */
struct Greedy final : OldShapeAlgorithm {
  Temple::JSF64 prng;

  double shape(const PositionCollection& positions, Shape shape) final {
    const unsigned N = positions.cols();
    const unsigned steps = 4 * facs.at(N - 2);

    Energy energy {positions, shapeCoordinates(shape)};
    auto permutation = Temple::iota<Vertex>(positions.cols());
    double lowestEnergy = std::numeric_limits<double>::max();

    bool foundBetterPermutation;
    for(unsigned m = 0; m < steps; ++m) {
      Temple::Random::shuffle(permutation, prng);
      double minimizationEnergy = energy(permutation);

      do {
        foundBetterPermutation = false;

        // Single swaps
        for(unsigned i = 0; i < N - 1 && m < steps; ++i) {
          for(unsigned j = i + 1; j < N && m < steps; ++j) {
            std::swap(permutation.at(i), permutation.at(j));

            double value = energy(permutation);
            ++m;
            if(value < minimizationEnergy) {
              minimizationEnergy = value;
              foundBetterPermutation = true;
              break;
            }

            // Swap back
            std::swap(permutation.at(i), permutation.at(j));
          }
        }

        // Propose double-swaps
        for(unsigned i = 0; i < N && !foundBetterPermutation; ++i) {
          for(unsigned j = i + 1; j < N && !foundBetterPermutation; ++j) {
            for(unsigned k = 0; k < N && k != i && !foundBetterPermutation; ++k) {
              for(unsigned l = k + 1; l < N && l != j; ++l) {
                std::swap(permutation.at(i), permutation.at(j));
                std::swap(permutation.at(k), permutation.at(l));

                // Calculate the value
                double value = energy(permutation);
                if(value < minimizationEnergy) {
                  minimizationEnergy = value;
                  foundBetterPermutation = true;
                  break;
                }
                ++m;
                if(m >= steps) {
                  return std::min(lowestEnergy, minimizationEnergy);
                }

                // Swap back
                std::swap(permutation.at(k), permutation.at(l));
                std::swap(permutation.at(i), permutation.at(j));
              }
            }
          }
        }

      } while(foundBetterPermutation);

      lowestEnergy = std::min(lowestEnergy, minimizationEnergy);
    }

    return lowestEnergy;
  }

  std::string name() const final {
    return "Greedy";
  }
};

/**
 * @brief Fixed steps steepest descent minimizer with shuffles
 *
 * Finds best swap to reduce permutation the most at each position and
 * minimizes. Shuffles when it's at a minimum. Unreliable.
 *
 * Complexity: Needs at least Theta(N - 2)! steps to discover enough state
 * space that it might find the minimum.
 */
struct SteepestDescent final : OldShapeAlgorithm {
  Temple::JSF64 prng;

  double shape(const PositionCollection& positions, Shape shape) final {
    const unsigned N = positions.cols();
    const unsigned steps = 4 * facs.at(N - 2);

    Energy energy {positions, shapeCoordinates(shape)};
    auto bestCandidate = Temple::iota<Vertex>(positions.cols());
    auto permutation = bestCandidate;
    double lowestEnergy = std::numeric_limits<double>::max();

    bool foundBetterPermutation;
    for(unsigned m = 0; m < steps; ++m) {
      Temple::Random::shuffle(bestCandidate, prng);
      double minimizationEnergy = energy(bestCandidate);
      lowestEnergy = std::min(lowestEnergy, minimizationEnergy);

      do {
        permutation = bestCandidate;
        foundBetterPermutation = false;

        // Propose double-swaps
        for(unsigned i = 0; i < N && !foundBetterPermutation; ++i) {
          for(unsigned j = i + 1; j < N && !foundBetterPermutation; ++j) {
            for(unsigned k = 0; k < N && !foundBetterPermutation; ++k) {
              for(unsigned l = k + 1; l < N; ++l) {
                std::swap(permutation.at(i), permutation.at(j));
                std::swap(permutation.at(k), permutation.at(l));

                // Calculate the value
                double value = energy(permutation);
                ++m;
                if(value < minimizationEnergy) {
                  minimizationEnergy = value;
                  foundBetterPermutation = true;
                  bestCandidate = permutation;
                }
                if(m >= steps) {
                  return std::min(lowestEnergy, minimizationEnergy);
                }

                // Swap back
                std::swap(permutation.at(k), permutation.at(l));
                std::swap(permutation.at(i), permutation.at(j));
              }
            }
          }
        }
      } while(foundBetterPermutation);

      lowestEnergy = std::min(lowestEnergy, minimizationEnergy);
    }

    return lowestEnergy;
  }

  std::string name() const final {
    return "Steepest";
  }
};

/**
 * @brief Aligns four atoms, then greedily zips the rest
 *
 * For each tuple of four atoms, aligns the positions, then greedily chooses
 * the best next sequence alignment.
 *
 * Really good for small distortions (<= 0.5). Suffers from its greediness
 * afterwards. Some branching might do it good.
 *
 * Complexity: Theta(N! / (N - 5)!) quaternion fits
 *
 * Potentially there is another optimization possible that could reduce the
 * scaling but there is a tradeoff involved that might counteract the formal
 * complexity decrease: Principally, spatial rotations of index mappings from
 * the underlying shape mappings need not be repeated. However, the question is
 * how to store them efficiently.
 */
struct AlignFive final : ShapeAlgorithm {
  using M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using IncrementalPermutation = std::unordered_map<unsigned, unsigned>;

  static Eigen::Matrix3d fitQuaternion(
    const PositionCollection& stator,
    const PositionCollection& rotor,
    const IncrementalPermutation& p
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

  static Permutation flatten(const IncrementalPermutation& p) {
    const unsigned P = p.size();
    Permutation flat(P);
    for(unsigned i = 0; i < P; ++i) {
      flat.at(i) = p.at(i);
    }
    return flat;
  }

  ResultPair narrow(
    const PositionCollection& stator,
    const PositionCollection& rotor,
    const std::set<Permutation>& correctRotations,
    bool inCorrectBranch,
    std::unordered_map<unsigned, unsigned> permutation,
    std::vector<unsigned> freeLeftVertices,
    std::vector<unsigned> freeRightVertices
  ) {
    const unsigned N = stator.cols();

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

    // Maps freeLeftVertices onto freeRightVertices
    auto subPermutation = Temple::iota<unsigned>(V);

    auto isCorrectBranch = [&](const auto& s) {
      return Temple::any_of(
        correctRotations,
        [&](const auto& rot) -> bool {
          // Rotation has to match the current permutation
          if(!Temple::all_of(
            permutation,
            [&](const unsigned left, const unsigned right) -> bool {
              return rot.at(left) == right;
            }
          )) {
            return false;
          }

          // Rotation has to match the sub-permutation too
          for(unsigned i = 0; i < V; ++i) {
            unsigned left = freeLeftVertices.at(i);
            unsigned right = freeRightVertices.at(s.at(i));
            if(rot.at(left) != right) {
              return false;
            }
          }

          return true;
        }
      );
    };

    decltype(subPermutation) bestPermutation;
    double minimalCost = std::numeric_limits<double>::max();
    do {
      double cost = 0.0;
      for(unsigned i = 0; i < V; ++i) {
        cost += costs(i, subPermutation.at(i));
      }

      if(cost < minimalCost) {
        if(inCorrectBranch && !bestPermutation.empty() && isCorrectBranch(bestPermutation) && !isCorrectBranch(subPermutation)) {
          std::cerr << "Excluding correct branch, decreasing cost from " << minimalCost << " to " << cost << " with non-optimum sub-permutation\n";
        }

        minimalCost = cost;
        bestPermutation = subPermutation;
      } else if(inCorrectBranch && isCorrectBranch(subPermutation)) {
        std::cerr << "Excluding correct branch due to cost " << cost << " > " << minimalCost << "\n";
      }
    } while(Temple::InPlace::next_permutation(subPermutation));

    // Fuse permutation and subpermutation
    for(unsigned i = 0; i < V; ++i) {
      permutation.emplace(freeLeftVertices.at(i), freeRightVertices.at(bestPermutation.at(i)));
    }

    assert(permutation.size() == N);

    auto R = fitQuaternion(stator, rotor, permutation);
    auto rotated = R * rotor;

    const double energy = Temple::accumulate(
      Temple::Adaptors::range(N),
      0.0,
      [&](const double carry, const unsigned i) -> double {
        return carry + (
          stator.col(i) - rotated.col(permutation.at(i))
        ).squaredNorm();
      }
    );

    return ResultPair {energy, flatten(permutation)};
  }

  ResultPair shape(
    const PositionCollection& positions,
    Shape shape,
    const Permutation& correctPermutation
  ) final {
    const unsigned N = positions.cols();
    if(N <= 4) {
      return AllPermutations {}.shape(positions, shape);
    }

    const auto shapeCoords = shapeCoordinates(shape);

    ResultPair minimal {std::numeric_limits<double>::max(), {}};
    std::unordered_map<unsigned, unsigned> permutation;

    const auto correctRotations = generateRotations(correctPermutation, shape);

    // i != j != k != l != m with {i, j, k, l, m} in [0, N)
    Temple::Loops::different(
      [&](const std::vector<unsigned>& indices) {
        permutation.clear();
        for(unsigned i = 0; i < 5; ++i) {
          permutation.emplace(i, indices[i]);
        }

        const unsigned i = indices[0];
        const unsigned j = indices[1];
        const unsigned k = indices[2];
        const unsigned l = indices[3];
        const unsigned m = indices[4];

        Eigen::Matrix3d R = fitQuaternion(positions, shapeCoords, permutation);
        auto rotatedShape = R * shapeCoords;

        double penalty = (
          (positions.col(0) - rotatedShape.col(i)).squaredNorm()
          + (positions.col(1) - rotatedShape.col(j)).squaredNorm()
          + (positions.col(2) - rotatedShape.col(k)).squaredNorm()
          + (positions.col(3) - rotatedShape.col(l)).squaredNorm()
          + (positions.col(4) - rotatedShape.col(m)).squaredNorm()
        );

        bool inCorrectBranch = Temple::any_of(
          correctRotations,
          [&](const auto& rot) {
            return rot[0] == i && rot[1] == j && rot[2] == k && rot[3] == l && rot[4] == m;
          }
        );
        if(inCorrectBranch) {
          std::cerr << "At correct permutation quintuple!\n";
        }

        if(penalty > minimal.first) {
          return;
        }

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

        auto narrowed = narrow(
          positions,
          rotatedShape,
          correctRotations,
          inCorrectBranch,
          permutation,
          std::move(freeLeftVertices),
          std::move(freeRightVertices)
        );

        if(narrowed.first < minimal.first) {
          minimal = narrowed;
        }
      },
      5,
      N
    );

    return minimal;
  }

  std::string name() const final {
    return "AlignFive";
  }
};

void writeEnergyStatistics() {
  /* Write an R file with all of the energy values for a particular shape
   */
  std::ofstream rFile("energy_statistics.R");
  std::vector<Shape> shapes;
  unsigned shapeNumber = 1;
  for(const Shape shape : allShapes) {
    if(size(shape) <= 4) {
      continue;
    }

    if(size(shape) > 8) {
      break;
    }

    std::cout << "Number of permutations for " << name(shape) << " shape: " << facs.at(size(shape)) << "\n";

    shapes.push_back(shape);

    auto shapeCoordinates = normalize(addOrigin(coordinates(shape)));
    const unsigned N = shapeCoordinates.cols();

    auto distorted = shapeCoordinates;
    const double distortionNorm = 0.2;
    distort(distorted, distortionNorm);

    Energy energy {shapeCoordinates, distorted};

    rFile << "shape" << shapeNumber << " <- c(";
    auto permutation = Temple::iota<Vertex>(N);
    rFile << energy(permutation);
    while(Temple::InPlace::next_permutation(permutation)) {
      rFile << ", " << energy(permutation);
    }
    rFile << ")\n";

    ++shapeNumber;
  }

  rFile << "shapeNames <- c(" << Temple::condense(Temple::map(shapes, [](auto s) { return "\"" + name(s) + "\""; })) << ")\n";
}

template<typename PRNG>
PositionCollection shuffle(const PositionCollection& positions, PRNG& prng) {
  const unsigned C = positions.cols();
  auto permutation = Temple::iota<Vertex>(C);
  Temple::Random::shuffle(permutation, prng);
  PositionCollection shuffled(3, C);
  for(unsigned i = 0; i < C; ++i) {
    shuffled.col(permutation.at(i)) = positions.col(i);
  }
  return shuffled;
}

namespace color {

std::ostream& boldMagenta(std::ostream& os) {
  os << "\033[1;35m";
  return os;
}

std::ostream& red(std::ostream& os) {
  os << "\033[31m";
  return os;
}

std::ostream& green(std::ostream& os) {
  os << "\033[32m";
  return os;
}

std::ostream& reset(std::ostream& os) {
  os << "\033[0m";
  return os;
}

} // namespace color

int main(int argc, char* argv[]) {
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help,h", "Produce help message")
    (
      "prng,p",
      boost::program_options::value<int>(),
      "Seed to initialize PRNG with."
    )
    (
      "shape,s",
      boost::program_options::value<unsigned>(),
      "Shape to run algorithms with"
    )
    (
      "distortion,d",
      boost::program_options::value<double>(),
      "Distortion norm to apply to shapes."
    )
    (
      "repeats,r",
      boost::program_options::value<unsigned>(),
      "Number of repetitions"
    )
  ;

  /* Parse */
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::command_line_parser(argc, argv).
    options(options_description).
    style(
      boost::program_options::command_line_style::unix_style
      | boost::program_options::command_line_style::allow_long_disguise
    ).run(),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);

  if(options_variables_map.count("help") > 0) {
    std::cout << options_description << "\n";
    return 0;
  }

  using namespace std::chrono;

  Temple::JSF64 prng;
  if(options_variables_map.count("seed")) {
    const int seed = options_variables_map["seed"].as<int>();
    prng.seed(seed);
    std::cout << "PRNG seeded from parameters: " << seed << ".\n";
  } else {
    std::random_device randomDevice;
    const int seed = std::random_device {}();
    std::cout << "PRNG seeded from random_device: " << seed << ".\n";
    prng.seed(seed);
  }

  std::vector<std::unique_ptr<ShapeAlgorithm>> algorithmPtrs;
  // algorithmPtrs.emplace_back(std::make_unique<Anneal>());
  // algorithmPtrs.emplace_back(std::make_unique<Tunnel>());
  // algorithmPtrs.emplace_back(std::make_unique<ThermodynamicAnneal>());
  // algorithmPtrs.emplace_back(std::make_unique<Greedy>());
  // algorithmPtrs.emplace_back(std::make_unique<SteepestDescent>());
  // algorithmPtrs.emplace_back(std::make_unique<CentroidElimination>());
  algorithmPtrs.emplace_back(std::make_unique<AlignFive>());

  const unsigned nameColWidth = 12;
  const unsigned timeColWidth = 6;
  const unsigned repeats = (options_variables_map.count("repeats") == 1)
    ? options_variables_map["repeats"].as<unsigned>()
    : 10;
  const double distortionNorm = (options_variables_map.count("distortion") == 1)
    ? options_variables_map["distortion"].as<double>()
    : 0.1;
  Shape shape = Shape::CappedSquareAntiprism;
  if(options_variables_map.count("shape") == 1) {
    unsigned shapeIndex = options_variables_map["shape"].as<unsigned>();
    if(shapeIndex < nShapes) {
      shape = static_cast<Shape>(shapeIndex);
    }
  }

  std::cout << "For shape " << name(shape) << " (size " << size(shape) << ") with distortion norm " << distortionNorm << "\n\n";

  AllPermutations referenceAlgorithm;
  std::cout << std::setw(nameColWidth) << referenceAlgorithm.name() << std::setw(timeColWidth) << "msec";

  for(auto& algorithmPtr : algorithmPtrs) {
    std::cout << std::setw(nameColWidth) << algorithmPtr->name() << std::setw(timeColWidth) << "msec" << std::setw(timeColWidth) << "S";
  }
  std::cout << "\n";

  const unsigned A = algorithmPtrs.size();
  std::vector<
    std::vector<double>
  > errors (A);

  std::vector<
    std::vector<double>
  > latencies (A);

  std::vector<
    std::vector<double>
  > speedups (A);

  for(unsigned i = 0; i < repeats; ++i) {
    auto coordinates = ShapeAlgorithm::shapeCoordinates(shape);
    distort(coordinates, distortionNorm);
    coordinates = normalize(shuffle(coordinates, prng));

    time_point<steady_clock> start, end;
    start = steady_clock::now();
    const auto referencePair = referenceAlgorithm.shape(coordinates, shape);
    end = steady_clock::now();
    const unsigned referenceLatency = duration_cast<microseconds>(end - start).count();
    std::cout << std::setw(nameColWidth) << referencePair.first << std::setw(timeColWidth) << (referenceLatency / 1000);
    auto rotations = ShapeAlgorithm::generateRotations(referencePair.second, shape);

    for(unsigned j = 0; j < A; ++j) {
      auto& algorithmPtr = algorithmPtrs.at(j);
      start = steady_clock::now();
      const auto valuePair = algorithmPtr->shape(coordinates, shape, referencePair.second);
      end = steady_clock::now();

      if(valuePair.second != referencePair.second) {
        std::cerr << "\n" << std::setw(25) << "Reference permutation" << " " << Temple::condense(referencePair.second) << "\n"
          << std::setw(25) << algorithmPtr->name() << " " << Temple::condense(valuePair.second) << "\n";
        if(rotations.count(valuePair.second) == 1) {
          std::cerr << "Algorithm result is a shape rotation of the reference permutation\n";
        } else {
          std::cerr << "Algorithm result is not a shape rotation of the reference permutation:\n";
          for(const auto& rotation : rotations) {
            std::cerr << "Rotation " << Temple::condense(rotation) << "\n";
          }
        }
      }

      if(valuePair.first < referencePair.first - 1e-5) {
        std::cout << color::boldMagenta << std::setw(nameColWidth) << valuePair.first << color::reset;
      } else if(std::fabs(valuePair.first - referencePair.first) < 1e-5) {
        std::cout << color::green << std::setw(nameColWidth) << valuePair.first << color::reset;
      } else {
        std::cout << color::red << std::setw(nameColWidth) << valuePair.first << color::reset;
      }
      errors.at(j).push_back(std::fabs(referencePair.first - valuePair.first));

      unsigned latency = duration_cast<microseconds>(end - start).count();
      std::cout << std::setw(timeColWidth) << (latency / 1000);
      latencies.at(j).push_back(latency);

      if(latency == 0) {
        std::cout << std::setw(timeColWidth) << "inf";
      } else {
        std::cout << std::setw(timeColWidth) << std::round(static_cast<double>(referenceLatency) / latency);
      }
      speedups.at(j).push_back(static_cast<double>(referenceLatency) / latency);
    }
    std::cout << "\n";
  }

  std::cout << "\n";
  std::cout << std::setw(nameColWidth + timeColWidth) << "averages";
  for(unsigned j = 0; j < A; ++j) {
    const double errorAverage = Temple::average(errors.at(j));
    if(errorAverage < 1e-10) {
      std::cout << color::green << std::setw(nameColWidth) << 0 << color::reset;
    } else {
      std::cout << color::red << std::setw(nameColWidth) << errorAverage << color::reset;
    }

    std::cout << std::setw(timeColWidth) << static_cast<unsigned>(Temple::average(latencies.at(j)) / 1000);
    std::cout << std::setw(timeColWidth) << static_cast<unsigned>(Temple::average(speedups.at(j)));
  }
  std::cout << "\n";
}
