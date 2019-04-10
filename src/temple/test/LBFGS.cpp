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

BOOST_AUTO_TEST_CASE(LBFGSTests) {
  auto const gradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::Ref<Eigen::VectorXd> gradients) {
    const double firstBracket = parameters[0] + 2 * parameters[1] - 7;
    const double secondBracket = 2 * parameters[0] + parameters[1] - 5;

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

  unsigned cycles = optimizer.optimize(positions, gradientTestFunction, gradientChecker);

  BOOST_CHECK(cycles < 100);
  BOOST_CHECK(std::fabs(positions[0] - 1.0) < 1e-3);
  BOOST_CHECK(std::fabs(positions[1] - 3.0) < 1e-3);
}
