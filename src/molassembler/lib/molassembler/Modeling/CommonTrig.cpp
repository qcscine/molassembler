// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Modeling/CommonTrig.h"

#include "dlib/optimization.h"

namespace molassembler {

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

namespace detail {

using DlibVector = ::dlib::matrix<double, 0, 1>;

struct SqrtLessDihedralLength {
  double operator () (const DlibVector& parameters) const {
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

    using namespace std;
    return (
      pow(a * sin(alpha) - c * cos(phi) * sin(beta), 2)
      + pow(c * sin(beta) * sin(phi), 2)
      + pow(a * cos(alpha) - b + c * cos(beta), 2)
    );
  }
};

struct SqrtLessDihedralLengthDerivative {
  DlibVector operator () (const DlibVector& parameters) const {
    const double& a = parameters(0);
    const double& b = parameters(1);
    const double& c = parameters(2);
    const double& alpha = parameters(3);
    const double& beta = parameters(4);
    const double& phi = parameters(5);

    using namespace std;
    // Figure out the derivatives per parameter
    DlibVector derivative(6);

    derivative(0) = 2 * (a + cos(alpha) * (-b + c * cos(beta)) -c * cos(phi) * sin(alpha) * sin(beta));
    derivative(1) = 2 * (b - a *cos(alpha) - c * cos(beta));
    derivative(2) = 2 * (c + (-b + a * cos(alpha)) * cos(beta) - a * cos(phi) * sin(alpha) * sin(beta));
    derivative(3) = 2 * a * (b - c * cos(beta)) * sin(alpha) - 2 * a * c * cos(alpha) * cos(phi) * sin(beta);
    derivative(4) = -2 * c * (a * cos(beta) * cos(phi) * sin(alpha) + (-b + a * cos(alpha)) * sin(beta));
    derivative(5) = 2 * a * c * sin(alpha) * sin(beta) * sin(phi);

    return derivative;
  }
};

} // namespace detail

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
   * For the limited value bounds that are commonly evaluated in Distance Geometry,
   * guessing cases can work, but generatilty is more easily achieved through numerics.
   *
   * Guessing combinations is slighly complex when it comes to the dihedral,
   * since this variable's bounds are not clamped to [0, pi] as would make
   * evaluating the min/max here easier.
   *
   * The dependence on the dihedral value is
   *
   * sqrt(
   *   ...
   *   + (factors > 0) * (
   *     cos(alpha) * cos(beta)
   *     - sin(alpha) * sin(beta) * cos(dihedral)
   *   )
   * )
   *
   * and so the dihedral length is proportional to -cos(dihedral).
   */

  /*auto cosBounds = [](const ValueBounds angleBounds) -> ValueBounds {
    using boostInterval = boost::numeric::interval<
      double,
      boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
          boost::numeric::interval_lib::rounded_transc_std<double>
        >,
        boost::numeric::interval_lib::checking_base<double>
      >
    >;

    auto transformedInterval = boost::numeric::cos(
      boostInterval(angleBounds.lower, angleBounds.upper)
    );

    return {
      boost::numeric::lower(transformedInterval),
      boost::numeric::upper(transformedInterval)
    };
  };

  ValueBounds dihedralCosBounds = cosBounds(dihedralBounds);*/

  /* Since we subtract either dihedralCosBounds' value, we invert the
   * lower/upper usage for the determination of the overall ValueBounds.
   *
   * If an angle's bounds do not include M_PI / 2, the situation is slightly
   * more difficult.
   *
   * When searching for the lower bound on a side with angle bounds whose upper
   * bound is smaller than M_PI / 2, the upper bound on its length may generate
   * a shorter dihedral length.
   */

  /*return {
    dihedralLength(
      (alphaBounds.lower < M_PI / 2 ? aBounds.upper : aBounds.lower),
      bBounds.lower,
      (betaBounds.lower < M_PI / 2 ? cBounds.upper : cBounds.lower),
      alphaBounds.lower,
      betaBounds.lower,
      dihedralCosBounds.upper
    ),
    dihedralLength(
      (alphaBounds.upper < M_PI / 2 ? aBounds.lower : aBounds.upper),
      bBounds.upper,
      (alphaBounds.upper < M_PI / 2 ? cBounds.lower : cBounds.upper),
      alphaBounds.upper,
      betaBounds.upper,
      dihedralCosBounds.lower
    )
  };*/

  using DlibVector = dlib::matrix<double, 0, 1>;

  DlibVector lower(6);
  lower = aBounds.lower, bBounds.lower, cBounds.lower, alphaBounds.lower,
        betaBounds.lower, dihedralBounds.lower;
  DlibVector upper(6);
  upper = aBounds.upper, bBounds.upper, cBounds.upper, alphaBounds.upper,
        betaBounds.upper, dihedralBounds.upper;

  /* Start minimum searches in two positions for each:
   * minimum: at lower bounds and at the medians
   * maximum: at the medians and at the upper bounds
   */

  DlibVector start = (lower + upper) / 2;

  double maxFromMedian = dlib::find_max_box_constrained(
    dlib::bfgs_search_strategy(),
    dlib::objective_delta_stop_strategy(),
    detail::SqrtLessDihedralLength (),
    detail::SqrtLessDihedralLengthDerivative (),
    start,
    lower,
    upper
  );

  start = upper;

  double maxFromUpper = dlib::find_max_box_constrained(
    dlib::bfgs_search_strategy(),
    dlib::objective_delta_stop_strategy(),
    detail::SqrtLessDihedralLength (),
    detail::SqrtLessDihedralLengthDerivative (),
    start,
    lower,
    upper
  );

  start = (lower + upper) / 2;

  double minFromMedian = dlib::find_min_box_constrained(
    dlib::bfgs_search_strategy(),
    dlib::objective_delta_stop_strategy(),
    detail::SqrtLessDihedralLength (),
    detail::SqrtLessDihedralLengthDerivative (),
    start,
    lower,
    upper
  );

  start = lower;

  double minFromLower = dlib::find_min_box_constrained(
    dlib::bfgs_search_strategy(),
    dlib::objective_delta_stop_strategy(),
    detail::SqrtLessDihedralLength (),
    detail::SqrtLessDihedralLengthDerivative (),
    start,
    lower,
    upper
  );

  return {
    std::sqrt(std::min(minFromMedian, minFromLower)),
    std::sqrt(std::max(maxFromMedian, maxFromUpper))
  };
}

} // namespace CommonTrig

} // namespace molassembler
