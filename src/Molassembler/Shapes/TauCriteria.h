/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Tau criteria for differentiating symmetries from coordinates
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_TAU_CRITERIA_H
#define INCLUDE_MOLASSEMBLER_SHAPES_TAU_CRITERIA_H

#include "Molassembler/Temple/constexpr/Math.h"
#include "Molassembler/Export.h"
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <vector>

namespace Scine {
namespace Molassembler {
namespace Shapes {
namespace Detail {

constexpr unsigned binomial(const unsigned n, const unsigned k) {
  return Temple::Math::factorial(n) / (
    Temple::Math::factorial(k) * Temple::Math::factorial(n - k)
  );
}

double tauFourPrime(const std::vector<double>& angles) {
  assert(std::is_sorted(std::begin(angles), std::end(angles)));

  const double beta = angles.back();
  const double alpha = angles.at(angles.size() - 2);

  constexpr double theta = Temple::Math::toRadians(109.5);

  return (
    (beta - alpha) / (2 * M_PI - theta)
    + (M_PI - beta) / (M_PI - theta)
  );
}

double tauFive(const std::vector<double>& angles) {
  assert(std::is_sorted(std::begin(angles), std::end(angles)));

  return (
    angles.back() - angles.at(angles.size() - 2)
  ) / Temple::Math::toRadians(60.0);
}

} // namespace Detail

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
MASM_EXPORT double tau(const std::vector<double>& angles) {
  constexpr unsigned anglesInSymmetryOfSizeFour = Detail::binomial(4, 2);
  constexpr unsigned anglesInSymmetryOfSizeFive = Detail::binomial(5, 2);

  if(angles.size() == anglesInSymmetryOfSizeFour) {
    return Detail::tauFourPrime(angles);
  }

  if(angles.size() == anglesInSymmetryOfSizeFive) {
    return Detail::tauFive(angles);
  }

  throw std::invalid_argument(
    "This function can only calculate tau values for centrosymmetries of sizes "
    "four or five. Your supplied vector of angles does not have matching size"
  );
}

} // namespace Shapes
} // namespace Molassembler
} // namespace Scine

#endif
