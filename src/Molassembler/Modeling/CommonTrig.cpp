/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Modeling/CommonTrig.h"

#include "Molassembler/Temple/Optimization/Lbfgs.h"

namespace Scine {
namespace Molassembler {
namespace CommonTrig {

double dihedralLength(
  const double a,
  const double b,
  const double c,
  const double alpha,
  const double beta,
  const double dihedral
) {
  return sqrt(
    a * a
    + b * b
    + c * c
    + 2 * (
      - a * b * std::cos(alpha)
      - b * c * std::cos(beta)
      + a * c * (
        std::cos(alpha) * std::cos(beta)
        - std::sin(alpha) * std::sin(beta) * std::cos(dihedral)
      )
    )
  );
}

namespace {

void dihedralLengthFn(const Eigen::VectorXd& parameters, double& value, Eigen::Ref<Eigen::VectorXd> gradient) {
  /* 0 -> a (i-j)
   * 1 -> b (j-k)
   * 2 -> c (k-l)
   * 3 -> alpha (angle between a and b / i-j-k)
   * 4 -> beta (angle between b and c / j-k-l)
   * 5 -> phi (dihedral on i-j-k-l)
   */
  const double& a = parameters(0);
  const double& b = parameters(1);
  const double& c = parameters(2);
  const double& alpha = parameters(3);
  const double& beta = parameters(4);
  const double& phi = parameters(5);

  value = (
    std::pow(a * std::sin(alpha) - c * std::cos(phi) * std::sin(beta), 2)
    + std::pow(c * std::sin(beta) * std::sin(phi), 2)
    + std::pow(a * std::cos(alpha) - b + c * std::cos(beta), 2)
  );

  gradient(0) = 2 * (a + std::cos(alpha) * (-b + c * std::cos(beta)) -c * std::cos(phi) * std::sin(alpha) * std::sin(beta));
  gradient(1) = 2 * (b - a *std::cos(alpha) - c * std::cos(beta));
  gradient(2) = 2 * (c + (-b + a * std::cos(alpha)) * std::cos(beta) - a * std::cos(phi) * std::sin(alpha) * std::sin(beta));
  gradient(3) = 2 * a * (b - c * std::cos(beta)) * std::sin(alpha) - 2 * a * c * std::cos(alpha) * std::cos(phi) * std::sin(beta);
  gradient(4) = -2 * c * (a * std::cos(beta) * std::cos(phi) * std::sin(alpha) + (-b + a * std::cos(alpha)) * std::sin(beta));
  gradient(5) = 2 * a * c * std::sin(alpha) * std::sin(beta) * std::sin(phi);
}

struct GradientBasedChecker {
  unsigned iterLimit = 100;
  double gradientLimit = 1e-5;

  template<typename StepValues>
  bool shouldContinue(unsigned iteration, const StepValues& step) {
    return (
      iteration < iterLimit
      && step.gradients.current.norm() > gradientLimit
    );
  }
};

} // namespace

ValueBounds dihedralLengthBounds(
  const ValueBounds& aBounds,
  const ValueBounds& bBounds,
  const ValueBounds& cBounds,
  const ValueBounds& alphaBounds,
  const ValueBounds& betaBounds,
  const ValueBounds& dihedralBounds
) {
  /* The dihedral length is a fairly complex 6-dimensional function.
   *
   * For the limited value bounds that are commonly evaluated in Distance
   * Geometry, guessing cases can work, but generality is more easily achieved
   * through numerical optimization.
   */

  Temple::Lbfgs<> optimizer;
  using VectorType = typename Temple::Lbfgs<>::VectorType;


  VectorType lower(6);
  lower << aBounds.lower, bBounds.lower, cBounds.lower, alphaBounds.lower,
        betaBounds.lower, dihedralBounds.lower;
  VectorType upper(6);
  upper << aBounds.upper, bBounds.upper, cBounds.upper, alphaBounds.upper,
        betaBounds.upper, dihedralBounds.upper;

  const Temple::Lbfgs<>::Box box { lower, upper };

  /* Start minimum searches in two positions for each:
   * minimum: at lower bounds and at the medians
   * maximum: at the medians and at the upper bounds
   */

  VectorType parameters = (box.minima + box.maxima) / 2;
  auto result = optimizer.maximize(
    parameters,
    box,
    &dihedralLengthFn,
    GradientBasedChecker {}
  );
  const double maxFromMedian = result.value;

  parameters = box.maxima;
  result = optimizer.maximize(
    parameters,
    box,
    &dihedralLengthFn,
    GradientBasedChecker {}
  );
  const double maxFromUpper = result.value;

  parameters = (box.minima + box.maxima) / 2;
  result = optimizer.minimize(
    parameters,
    box,
    &dihedralLengthFn,
    GradientBasedChecker {}
  );
  const double minFromMedian = result.value;

  parameters = lower;
  result = optimizer.minimize(
    parameters,
    box,
    &dihedralLengthFn,
    GradientBasedChecker {}
  );
  const double minFromLower = result.value;

  return {
    std::sqrt(std::min(minFromMedian, minFromLower)),
    std::sqrt(std::max(maxFromMedian, maxFromUpper))
  };
}

} // namespace CommonTrig
} // namespace Molassembler
} // namespace Scine
