/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Tau criteria for differentiating symmetries from coordinates
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_TAU_CRITERIA_H
#define INCLUDE_MOLASSEMBLER_SHAPES_TAU_CRITERIA_H

#include "temple/constexpr/Math.h"
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <vector>

namespace Scine {

namespace Symmetry {

namespace detail {

constexpr unsigned binomial(const unsigned n, const unsigned k) {
  return temple::Math::factorial(n) / (
    temple::Math::factorial(k) * temple::Math::factorial(n - k)
  );
}

double tauFourPrime(const std::vector<double>& angles) {
  assert(std::is_sorted(std::begin(angles), std::end(angles)));

  const double beta = angles.back();
  const double alpha = angles.at(angles.size() - 2);

  constexpr double theta = temple::Math::toRadians(109.5);

  return (
    (beta - alpha) / (2 * M_PI - theta)
    + (M_PI - beta) / (M_PI - theta)
  );
}

double tauFive(const std::vector<double>& angles) {
  assert(std::is_sorted(std::begin(angles), std::end(angles)));

  return (
    angles.back() - angles.at(angles.size() - 2)
  ) / temple::Math::toRadians(60.0);
}

} // namespace detail

/**
 * @brief Calculates the tau value for four and five angle symmetries
 *
 * @complexity{@math{\Theta(1)}}
 *
 * @param angles A sorted vector of the angles in your central symmetry. This
 *   must be of size four or five.
 *
 * @return τ₄' or τ₅
 */
double tau(const std::vector<double>& angles) {
  constexpr unsigned anglesInSymmetryOfSizeFour = detail::binomial(4, 2);
  constexpr unsigned anglesInSymmetryOfSizeFive = detail::binomial(5, 2);

  if(angles.size() == anglesInSymmetryOfSizeFour) {
    return detail::tauFourPrime(angles);
  }

  if(angles.size() == anglesInSymmetryOfSizeFive) {
    return detail::tauFive(angles);
  }

  throw std::invalid_argument(
    "This function can only calculate tau values for centrosymmetries of sizes "
    "four or five. Your supplied vector of angles does not have matching size"
  );
}

} // namespace Symmetry

} // namespace Scine

#endif
