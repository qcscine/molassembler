/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Optimization/Lbfgs.h"

using namespace Scine::Molassembler;

template<typename FloatType>
struct GradientBasedChecker {
  unsigned iterLimit = 100;
  FloatType gradientLimit = 1e-5;

  template<typename StepValues>
  bool shouldContinue(unsigned iteration, const StepValues& step) {
    return (
      iteration < iterLimit
      && step.gradients.current.norm() > gradientLimit
    );
  }
};

BOOST_AUTO_TEST_CASE(LBFGSSimpleMinimization) {
  /* Booth function
   * Ellipsoid f(x, y) = (x + 2y - 7)² + (2x + y - 5)²
   * Minimum at x = 1, y = 3
   */
  const auto gradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::Ref<Eigen::VectorXd> gradients) {
    const double& x = parameters[0];
    const double& y = parameters[1];
    const double firstBracket = x + 2 * y - 7;
    const double secondBracket = 2 * x + y - 5;

    value = (
      std::pow(firstBracket, 2)
      + std::pow(secondBracket, 2)
    );
    gradients[0] = 2 * firstBracket + 4 * secondBracket;
    gradients[1] = 4 * firstBracket + 2 * secondBracket;
  };

  Temple::Lbfgs<double, 16> optimizer;
  GradientBasedChecker<double> gradientChecker;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.75 * M_PI;

  const auto result = optimizer.minimize(positions, gradientTestFunction, gradientChecker);

  BOOST_CHECK_MESSAGE(
    result.iterations < 100,
    "Expected convergence in less than 100 cycles, got " << result.iterations
  );
  BOOST_CHECK(std::fabs(positions[0] - 1.0) < 1e-3);
  BOOST_CHECK(std::fabs(positions[1] - 3.0) < 1e-3);
}

BOOST_AUTO_TEST_CASE(LBFGSSimpleMaximization) {
  /* Very simple parabola -((x-4)² + (y-2)²) + 4
   * Maximum at 4, 2
   */
  const auto gradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::Ref<Eigen::VectorXd> gradients) {
    const double& x = parameters[0];
    const double& y = parameters[1];

    value = -(std::pow(x - 4, 2) + std::pow(y - 2, 2)) + 4;
    gradients[0] = -2 * (x - 4);
    gradients[1] = -2 * (y - 2);
  };

  Temple::Lbfgs<double, 16> optimizer;
  GradientBasedChecker<double> gradientChecker;

  Eigen::VectorXd positions(2);
  positions << 2.0, -1.0;

  Eigen::VectorXd expectedMaximum(2);
  expectedMaximum << 4.0, 2.0;

  const auto result = optimizer.maximize(positions, gradientTestFunction, gradientChecker);

  BOOST_CHECK_MESSAGE(
    result.iterations < 100,
    "Expected convergence in less than 100 cycles, got " << result.iterations
  );

  for(unsigned i = 0; i < 2; ++i) {
    BOOST_CHECK_MESSAGE(
      std::fabs(positions[i] - expectedMaximum[i]) < 1e-3,
      "Position parameter " << i << " is not at the expected maximum."
      << " Expected " << expectedMaximum[i] << ", got " << positions[i]
    );
  }
}

BOOST_AUTO_TEST_CASE(LBFGSBraninMinimization) {
  /* f(x, y) = (-1.275 * x² / pi² + 4 x / pi + y - 6)² + (10 - 5 / (4 pi)) * cos(x) + 10
   * Minima at (2n pi, 2n pi)
   */
  const auto gradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::Ref<Eigen::VectorXd> gradients) {
    const double& x = parameters[0];
    const double& y = parameters[1];

    constexpr double alpha = -1.275 / (M_PI * M_PI);
    constexpr double beta = 4 / M_PI;
    constexpr double cosinePrefactor = (10 - 5 / (4 * M_PI));

    const double polynomialBracket = (alpha * x * x + beta * x + y - 6);

    value = (
      std::pow(polynomialBracket, 2)
      + cosinePrefactor * std::cos(x)
      + 10
    );

    gradients[0] = 2 * polynomialBracket * (2 * alpha * x + beta) - cosinePrefactor * std::sin(x);
    gradients[1] = 2 * polynomialBracket;
  };

  using OptimizerType = Temple::Lbfgs<double, 16>;
  OptimizerType optimizer;
  GradientBasedChecker<double> gradientChecker;
  Eigen::VectorXd positions(2);
  positions << M_PI - 0.1, M_PI - 0.1;

  constexpr double expectedResult = 0.39788735772973816;

  const auto result = optimizer.minimize(
    positions,
    gradientTestFunction,
    gradientChecker
  );

  BOOST_CHECK_MESSAGE(
    result.iterations < 100,
    "Expected convergence in less than 100 cycles, got " << result.iterations
  );

  BOOST_CHECK_MESSAGE(
    std::fabs(result.value - expectedResult) < 1e-3,
    "Expected minimum of value " << expectedResult << ", got " << result.value << " instead"
  );
}

BOOST_AUTO_TEST_CASE(LBFGSBoxedMinimization) {
  /* f(x, y) = - cos x - 0.5 cos y
   * Minima at (2n pi, 2n pi)
   * In box min [0.1, 0], max [pi, pi] minimum should be [0.1, 0]
   */
  const auto gradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::Ref<Eigen::VectorXd> gradients) {
    const double& x = parameters[0];
    const double& y = parameters[1];

    value = (
      - std::cos(x)
      - 0.5 * std::cos(y)
    );
    gradients[0] = std::sin(x);
    gradients[1] = 0.5 * std::sin(y);
  };

  using OptimizerType = Temple::Lbfgs<double, 16>;

  OptimizerType optimizer;
  GradientBasedChecker<double> gradientChecker;
  Eigen::VectorXd positions(2), boxMinima(2), boxMaxima(2);
  boxMinima << 0.1, 0;
  positions << 0.6, 0.5;
  boxMaxima << M_PI, M_PI;

  const typename OptimizerType::Box box {
    boxMinima,
    boxMaxima
  };

  const auto result = optimizer.minimize(
    positions,
    box,
    gradientTestFunction,
    gradientChecker
  );

  BOOST_CHECK_MESSAGE(
    result.iterations < 100,
    "Expected convergence in less than 100 cycles, got " << result.iterations
  );
  BOOST_CHECK_MESSAGE(
    std::fabs(positions[0] - boxMinima[0]) < 1e-3,
    "Expected x_min = x_box_min (0.1), but is " << positions[0] << " instead"
  );
  BOOST_CHECK_MESSAGE(
    std::fabs(positions[1] - boxMinima[1]) < 1e-3,
    "Expected y_min = y_box_min (0), but is " << positions[1] << " instead"
  );
}
