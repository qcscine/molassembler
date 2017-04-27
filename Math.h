#ifndef INCLUDE_CONSTEXPR_MAGIC_MATH_H
#define INCLUDE_CONSTEXPR_MAGIC_MATH_H

/* TODO 
 * -remove dependency
 */
#include "static_math/cmath.h"

namespace ConstexprMagic {

namespace Math {

// Implements the 
constexpr double asin (const double& x) {
  if(!(-1 < x && x < 1)) {
    throw "inverse sine domain error: only real if -1 < x < 1";
  }

  const double epsilon = 1e-8;
  double value = 0;
  double change = 1;

  double upper_factorial = 1;
  double lower_factorial = 1;

  for(unsigned n = 0; change > epsilon; n++) {
    if(n > 0) {
      upper_factorial *= 2 * (n - 1) + 1;
      lower_factorial *= 2 * n;
    }

    change = value;

    value += (
      upper_factorial / lower_factorial
    ) * smath::pow(x, 2 * n + 1);

    change = smath::abs(change - value);
  }

  return value;
}

constexpr double acos (const double& x) {
  if(!(-1 < x && x < 1)) {
    throw "inverse cosine domain error: only real if -1 < x < 1";
  }

  return M_PI / 2 - asin(x);
}

} // namespace Math

} // namespace ConstexprMagic

#endif
