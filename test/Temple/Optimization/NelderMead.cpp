/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Optimization/NelderMead.h"

using namespace Scine::Molassembler;

struct NelderMeadHimmelblau {
  double operator() (const Eigen::VectorXd& parameters) {
    assert(parameters.size() == 2);

    const double x = parameters(0);
    const double y = parameters(1);

    const double firstBracket = (x * x + y - 11);
    const double secondBracket = (x + y * y - 7);

    return firstBracket * firstBracket + secondBracket * secondBracket;
  }
};

BOOST_AUTO_TEST_CASE(NelderMead, *boost::unit_test::label("Temple")) {
  Eigen::Matrix<double, 2, 3> simplexVertices;
  simplexVertices << -3.0,  0.0,  0.0,
                     -1.5,  0.0, -3.0;

  struct NelderMeadChecker {
    bool shouldContinue(unsigned iteration, double /* lowestValue */, double stddev) {
      return iteration < 1000 && stddev > 0.01;
    }
  };

  auto optimizationResult = Temple::NelderMead<>::minimize(
    simplexVertices,
    NelderMeadHimmelblau {},
    NelderMeadChecker {}
  );

  BOOST_CHECK_MESSAGE(
    std::fabs(optimizationResult.value) <= 0.1,
    "Nelder-Mead does not find minimization of Himmelblau function, value is "
    << optimizationResult.value << " after " << optimizationResult.iterations
    << " iterations. Simplex vertices are: " << simplexVertices
  );
}
