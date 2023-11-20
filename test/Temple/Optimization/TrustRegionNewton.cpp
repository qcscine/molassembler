/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Optimization/TrustRegion.h"
using namespace Scine::Molassembler;

struct Himmelblau {
  double operator() (const Eigen::VectorXd& parameters) {
    assert(parameters.size() == 2);

    const double x = parameters(0);
    const double y = parameters(1);

    const double firstBracket = (x * x + y - 11);
    const double secondBracket = (x + y * y - 7);

    return firstBracket * firstBracket + secondBracket * secondBracket;
  }

  void operator() (
    const Eigen::VectorXd& parameters,
    double& value,
    Eigen::Ref<Eigen::VectorXd> gradient,
    Eigen::Ref<Eigen::MatrixXd> hessian
  ) {
    assert(parameters.size() == 2);
    assert(gradient.size() == 2);
    assert(hessian.cols() == 2 && hessian.rows() == 2);

    const double x = parameters(0);
    const double y = parameters(1);

    const double firstBracket = (x * x + y - 11);
    const double secondBracket = (x + y * y - 7);

    value = firstBracket * firstBracket + secondBracket * secondBracket;

    gradient(0) = 4 * x * firstBracket + 2 * secondBracket;
    gradient(1) = 4 * y * secondBracket + 2 * firstBracket;

    hessian(0, 0) = 8 * x * x + 4 * firstBracket + 2;
    hessian(1, 1) = 8 * y * y + 4 * secondBracket + 2;
    hessian(0, 1) = 4 * x + 4 * y;
    hessian(1, 0) = hessian(0, 1);
  }

  bool shouldContinue(
    const unsigned iteration,
    const double /* value */,
    const Eigen::VectorXd& gradient
  ) {
    return (
      iteration <= 1000
      && gradient.squaredNorm() > 1e-3
    );
  }
};

BOOST_AUTO_TEST_CASE(TrustRegionNewton, *boost::unit_test::label("Temple")) {
  Eigen::VectorXd parameters = Eigen::VectorXd::Random(2);
  Eigen::VectorXd passParameters = parameters;

  auto optimizationResult = Temple::TrustRegionOptimizer<>::minimize(
    parameters,
    Himmelblau {},
    Himmelblau {}
  );

  BOOST_CHECK_MESSAGE(
    std::fabs(optimizationResult.value) <= 1e-5,
    "Newton-Raphson trust region does not find minimization of Himmelblau function, value is "
    << optimizationResult.value << " after " << optimizationResult.iterations
    << " iterations at " << parameters.transpose() << ". Gradient norm is " << optimizationResult.gradient.norm()
  );
}
