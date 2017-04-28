#ifndef INCLUDE_CONSTEXPR_MAGIC_MATH_H
#define INCLUDE_CONSTEXPR_MAGIC_MATH_H

#include <cmath>

/* TODO 
 * - Add periodicities of the trigonometric functions
 * - Investigate ill-conditioned quality of inverse trig functions at specific
 *   angles (see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions)
 * - Improve power function, it's probably hella awful
 */
namespace ConstexprMagic {

namespace Math {

// Not a whole lot you can do wrong here
template<typename T>
inline constexpr T abs(const T& x) {
  return (x >= 0) ? x : -x;
}

/* Integer version just calls the unsigned power function
 * TODO lots can go wrong here!
 */
constexpr double pow(const double& base, const int& exponent) {
  if(exponent < 0) {
    return 1 / pow(base, static_cast<unsigned>(abs(exponent)));
  }

  else return pow(base, static_cast<unsigned>(exponent));
}

// Really weak first implementation
constexpr double pow(const double& base, const unsigned& exponent) {
  double value = base;

  for(unsigned n = 1; n < exponent; n++) {
    value *= base;
  }

  return value;
}

template<typename T>
inline constexpr auto floor(const T& x) -> decltype(std::floor(x)) {
  return (
    (int(x) == x) ? int(x) : (
      (x >= 0.0) ? int(x) : int(x) - 1
    )
  );
}

/* Implements Newton's iteration to compute the square root of a positive number
 */
constexpr double sqrt(const double& x) {
  if(x < 0) {
    throw "Square-root domain error: Only real if x >= 0!";
  }

  const double epsilon = 1e-8;
  double value = 1;
  double previous = 2;

  while(abs(previous - value) > epsilon) {
    // store the previous value
    previous = value;

    // compute next iteration
    value = 0.5 * (value + x / value);
  }

  return value;
}

/* Implements the infinite series where the derivative is expanded as a binomial
 * series and every term is integrated. Wikipedia says this implementation is 
 * probably ill-conditioned at some angles 
 */
constexpr double asin(const double& x) {
  if(!(-1 < x && x < 1)) {
    throw "Inverse sine domain error: only real if -1 < x < 1!";
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
    ) * pow(x, 2 * n + 1);

    change = abs(change - value);
  }

  return value;
}

constexpr double acos(const double& x) {
  if(!(-1 < x && x < 1)) {
    throw "Inverse cosine domain error: only real if -1 < x < 1!";
  }

  return M_PI / 2 - asin(x);
}

} // namespace Math

} // namespace ConstexprMagic

#endif
