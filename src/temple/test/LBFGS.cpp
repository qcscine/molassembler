/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "temple/LBFGS.h"

template<typename FloatType>
struct GradientBasedChecker {
  unsigned iterLimit = 100;
  FloatType gradientLimit = 1e-5;

  bool checkMaxIterations(unsigned iteration) {
    return iteration >= iterLimit;
  }

  bool checkConvergence(const Eigen::VectorXd& /* parameters */, const double /* value */, const Eigen::VectorXd& gradients) {
    return gradients.norm() <= gradientLimit;
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

  temple::LBFGS<double, 16> optimizer;
  GradientBasedChecker<double> gradientChecker;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.75 * M_PI;

  const unsigned cycles = optimizer.minimize(positions, gradientTestFunction, gradientChecker);

  BOOST_CHECK_MESSAGE(
    cycles < 100,
    "Expected convergence in less than 100 cycles, got " << cycles
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

  temple::LBFGS<double, 16> optimizer;
  GradientBasedChecker<double> gradientChecker;

  Eigen::VectorXd positions(2);
  positions << 2.0, -1.0;

  Eigen::VectorXd expectedMaximum(2);
  expectedMaximum << 4.0, 2.0;

  const unsigned cycles = optimizer.maximize(positions, gradientTestFunction, gradientChecker);

  BOOST_CHECK_MESSAGE(
    cycles < 100,
    "Expected convergence in less than 100 cycles, got " << cycles
  );

  for(unsigned i = 0; i < 2; ++i) {
    BOOST_CHECK_MESSAGE(
      std::fabs(positions[i] - expectedMaximum[i]) < 1e-3,
      "Position parameter " << i << " is not at the expected maximum."
      << " Expected " << expectedMaximum[i] << ", got " << positions[i]
    );
  }
}

BOOST_AUTO_TEST_CASE(LBFGSCosineMinimization) {
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

  using OptimizerType = temple::LBFGS<double, 16>;
  OptimizerType optimizer;
  GradientBasedChecker<double> gradientChecker;
  Eigen::VectorXd positions(2);
  positions << M_PI - 0.1, M_PI - 0.1;

  Eigen::VectorXd expectedMinimum(2);
  expectedMinimum << 0, 0;

  const unsigned cycles = optimizer.maximize(
    positions,
    gradientTestFunction,
    gradientChecker
  );

  Eigen::VectorXd tmpGrad(2);
  double value;
  gradientTestFunction(positions, value, tmpGrad);
  std::cout << "Minimized to " << positions.transpose() << " in " << cycles << " cycles. Value is " << value << "\n";

  BOOST_CHECK_MESSAGE(
    cycles < 100,
    "Expected convergence in less than 100 cycles, got " << cycles
  );

  for(unsigned i = 0; i < 2; ++i) {
    BOOST_CHECK_MESSAGE(
      std::fabs(positions[i] - expectedMinimum[i]) < 1e-3,
      "Position parameter " << i << " is not at the expected maximum."
      << " Expected " << expectedMinimum[i] << ", got " << positions[i]
    );
  }
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

  using OptimizerType = temple::LBFGS<double, 16>;

  OptimizerType optimizer;
  GradientBasedChecker<double> gradientChecker;
  Eigen::VectorXd positions(2), boxMinima(2), boxMaxima(2);
  boxMinima << 0.1, 0;
  positions << M_PI - 0.1, M_PI - 0.1;
  boxMaxima << M_PI, M_PI;

  const typename OptimizerType::Box box {
    boxMinima,
    boxMaxima
  };

  const unsigned cycles = optimizer.minimize(
    positions,
    box,
    gradientTestFunction,
    gradientChecker
  );

  Eigen::VectorXd tmpGrad(2);
  double value;
  gradientTestFunction(positions, value, tmpGrad);
  std::cout << "Box-minimized to " << positions.transpose() << " in " << cycles << " cycles. Value is " << value << "\n";

  BOOST_CHECK_MESSAGE(
    cycles < 100,
    "Expected convergence in less than 100 cycles, got " << cycles
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

BOOST_AUTO_TEST_CASE(LBFGSBoxedNonConstrainedMaximization) {
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

  using OptimizerType = temple::LBFGS<double, 16>;
  OptimizerType optimizer;
  GradientBasedChecker<double> gradientChecker;
  Eigen::VectorXd positions(2), boxMinima(2), boxMaxima(2);
  boxMinima << 0, 0;
  positions << 0, 4;
  boxMaxima << 6, 6;

  Eigen::VectorXd expectedMaximum(2);
  expectedMaximum << 4, 2;

  const typename OptimizerType::Box box {
    boxMinima,
    boxMaxima
  };

  const unsigned cycles = optimizer.maximize(
    positions,
    box,
    gradientTestFunction,
    gradientChecker
  );

  Eigen::VectorXd tmpGrad(2);
  double value;
  gradientTestFunction(positions, value, tmpGrad);
  std::cout << "Box-maximized to " << positions.transpose() << " in " << cycles << " cycles. Value is " << value << "\n";

  BOOST_CHECK_MESSAGE(
    cycles < 100,
    "Expected convergence in less than 100 cycles, got " << cycles
  );

  for(unsigned i = 0; i < 2; ++i) {
    BOOST_CHECK_MESSAGE(
      std::fabs(positions[i] - expectedMaximum[i]) < 1e-3,
      "Position parameter " << i << " is not at the expected maximum."
      << " Expected " << expectedMaximum[i] << ", got " << positions[i]
    );
  }
}

BOOST_AUTO_TEST_CASE(LBFGSBoxedConstrainingMaximization) {
  /* f(x, y) = - cos x - 0.5 cos y
   * Maxima at (n pi, n pi)
   * In box of min [0, 0], max [2 pi, pi - 0.1] max should be [pi, pi - 0.1]
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

  using OptimizerType = temple::LBFGS<double, 16>;
  OptimizerType optimizer;
  GradientBasedChecker<double> gradientChecker;
  Eigen::VectorXd positions(2), boxMinima(2), boxMaxima(2);
  boxMinima << 0, 0;
  positions << 0.1, 0.1;
  boxMaxima << 2 * M_PI, M_PI - 0.1;

  const typename OptimizerType::Box box {
    boxMinima,
    boxMaxima
  };

  const unsigned cycles = optimizer.maximize(
    positions,
    box,
    gradientTestFunction,
    gradientChecker
  );

  Eigen::VectorXd tmpGrad(2);
  double value;
  gradientTestFunction(positions, value, tmpGrad);
  std::cout << "Box-maximized to " << positions.transpose() << " in " << cycles << " cycles. Value is " << value << "\n";

  BOOST_CHECK_MESSAGE(
    cycles < 100,
    "Expected convergence in less than 100 cycles, got " << cycles
  );
  BOOST_CHECK_MESSAGE(
    std::fabs(positions[0] - M_PI) < 1e-3,
    "Expected x_min = pi, but is " << positions[0] << " instead"
  );
  BOOST_CHECK_MESSAGE(
    std::fabs(positions[1] - (M_PI - 0.1)) < 1e-3,
    "Expected y_min = y_max_boxed (pi - 0.1), but is " << positions[1] << " instead"
  );
}
