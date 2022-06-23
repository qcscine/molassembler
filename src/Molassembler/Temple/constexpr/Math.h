/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief constexpr math implementations
 *
 * Provides \c constexpr basic mathematical function implementations, some
 * logical functions and floating-point comparison helpers.
 *
 * @warning Do not use these math functions for anything other than constexpr
 * evaluations. The standard library implementations are most definitely
 * better. In many cases, even though it is not required by the standard, the
 * supplied STL functions are constexpr.
 */

#ifndef INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_MATH_H
#define INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_MATH_H

#include "Molassembler/Temple/Preprocessor.h"

#include <cmath>
#include <limits>
#include <type_traits>

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Math {
namespace Traits {

template<typename T, typename U>
using enableIfFloatingWithReturn = std::enable_if_t<
  std::is_floating_point<T>::value,
  U
>;

template<typename T, typename U>
using enableIfIntegralWithReturn = std::enable_if_t<
  std::is_integral<T>::value,
  U
>;

template<typename T, typename U>
using enableIfArithmeticWithReturn = std::enable_if_t<
  std::is_arithmetic<T>::value,
  U
>;

template<typename FloatingPoint>
PURITY_STRONG constexpr inline std::enable_if_t<
  std::is_floating_point<FloatingPoint>::value,
  bool
> isnan(const FloatingPoint x) {
  // see notes at http://en.cppreference.com/w/cpp/numeric/math/isnan
  return x != x;
}

} // namespace Traits

/* Logic */
//! Template parameter-pack exclusive or of booleans
template<typename ... Bools>
constexpr bool XOR(Bools ... bools);

/* Some very basic math functions for arithmetic types */
//! Absolute value
template<typename T>
inline constexpr Traits::enableIfArithmeticWithReturn<T, T> abs(T x) noexcept;

//! Maximum of two values
template<typename T>
constexpr Traits::enableIfArithmeticWithReturn<T, T> max(T a, T b) noexcept;

//! Minimum of two values
template<typename T>
constexpr Traits::enableIfArithmeticWithReturn<T, T> min(T a, T b) noexcept;

/* Floating-point math functions */

//! Convert angular degrees to radians
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> toRadians(T inDegrees) noexcept;

//! Convert angular radians to degrees
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> toDegrees(T inRadians) noexcept;

//! module function for arbitrary types
template<typename T, typename U>
constexpr Traits::enableIfFloatingWithReturn<T, T> fmod(T value, U divider) noexcept;

//! Ceiling function, no overflow check
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, long> ceil(T value) noexcept;

//! Floor function
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, int> floor(T value) noexcept;

//! Power of a number
template<typename T>
constexpr Traits::enableIfArithmeticWithReturn<T, T> pow(T base, unsigned exponent) noexcept;

template<typename T>
constexpr Traits::enableIfArithmeticWithReturn<T, double> pow(T base, int exponent) noexcept;

//! Square root
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> sqrt(T x);

//! Factorial
template<typename T>
constexpr Traits::enableIfIntegralWithReturn<T, T> factorial(T x);

//! Natural logarithm
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> ln(T x);

//! Base-10 logarithm
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> log10(T x);

//! Arbitrary base logarithm
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> log(T x, T base);

// Inverse trigonometry
/*! Computes the inverse sine function.
 *
 * NOTE: Accurate to only ~1e-9 absolute deviation close to domain boundaries
 */
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> asin(T x);

//! Inverse cosine
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> acos(T x);

//! Inverse tangens
template<typename T>
constexpr Traits::enableIfFloatingWithReturn<T, T> atan(T x);


/* Implementations begin here ------------------------------------------------*/

namespace Detail { // Implementation helpers

// Specialization of TPPSum for empty parameter pack
constexpr unsigned TPPSum() {
  return 0;
}

// Template parameter pack sum (needed for XOR function)
template<typename T1, typename... T>
constexpr unsigned TPPSum(T1 a, T ... pack) {
  return a + TPPSum(pack ...);
}

/* Based on series expansion of ln x:
 *
 *   y(x) = (x - 1) / (x + 1)
 *   ln x = 2 [ y + y^3/3 + y^5/5 + ...]
 *
 *   for Re x ≥ 0 and x ≠ 0
 */
template<typename T>
PURITY_STRONG constexpr T lnSeries(const T x) {
  if(x <= 0) {
    throw "Ln domain error: x <= 0";
  }

  const T epsilon = std::numeric_limits<T>::epsilon();

  // Since every term, y is multiplied by y², we keep a counting square
  T countingPower { // at n = 1
    (x - 1) / (x + 1)
  };

  // This is just y², which is what we multiply which each term
  const T iterationMultiplier = countingPower * countingPower;

  // Initial value for n = 1 (terms up until y)
  T value {
    2 * countingPower
  };

  // Term n = 0
  T previous {0};

  for(unsigned n = 3; Temple::Math::abs(previous - value) > epsilon; n += 2) {
    // previous iteration
    previous = value;

    // Update the power: y^(n-2) * y^2 becomes y^n
    countingPower *= iterationMultiplier;

    // Add the full term to the running value
    value += 2 * countingPower / n;
  }

  return value;
}

/* Implements an approximation to asin over 0 < x < 1 with accuracy > 2e-8
 * from Abramowitz, M., Stegun, I.: Handbook of Mathematical Functions, p. 64,
 * 1964, from http://people.math.sfu.ca/~cbm/aands/abramowitz_and_stegun.pdf
 */
template<typename T>
PURITY_STRONG constexpr T asinApprox(const T x) {
  if(!(0.0 <= x && x <= 1.0)) {
    throw "Asin approximation domain error: only applicable for 0 < x < 1!";
  }

  const T x2 = x * x;
  const T x4 = x2 * x2;

  return (
    M_PI / 2
    - Math::sqrt(1 - x) * (
      1.5707963050
      - 0.2145988016 * x
      + 0.0889789874 * x2
      - 0.0501743046 * x2 * x
      + 0.0308918810 * x4
      - 0.0170881256 * x4 * x
      + 0.0066700901 * x4 * x2
      - 0.0012624911 * x4 * x2 * x
    )
  );
}

} // namespace Detail

template<typename ... Bools>
constexpr bool XOR(Bools ... bools) {
  return Detail::TPPSum(bools ...) == 1;
}

template<typename T>
PURITY_STRONG inline constexpr Traits::enableIfArithmeticWithReturn<T, T> abs(const T x) noexcept {
  return (x >= 0) ? x : -x;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfArithmeticWithReturn<T, T> max(const T a, const T b) noexcept {
  return (a > b) ? a : b;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfArithmeticWithReturn<T, T> min(const T a, const T b) noexcept {
  return (a < b) ? a : b;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> toRadians(const T inDegrees) noexcept {
  return M_PI * inDegrees / 180;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> toDegrees(const T inRadians) noexcept {
  return 180 * inRadians / M_PI;
}

template<typename T, typename U>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T>  fmod(const T value, const U divider) noexcept {
 T remainder = value;
 while (divider < remainder) {
   remainder -= divider;
 }
 return remainder;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, long> ceil(const T value) noexcept {
  // Truncate to an integer
  const auto truncated = static_cast<long>(value);

  // we first check if the given floating point number is actually an integer 
  // this avoids rounding mistakes, occuring with some compilers
  const double eps = 1e-12;
  if (fmod(value, 1) < eps && abs(truncated - value) < eps) { 
    return truncated;
  }

  if(truncated < value) {
    return truncated + 1;
  }

  return truncated;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, int> floor(const T value) noexcept {
  // Truncate to an int
  const auto truncated = static_cast<int>(value);

  if(truncated > value) {
    return truncated - 1;
  }

  return truncated;
}

// Really weak first implementation
template<typename T>
PURITY_STRONG constexpr Traits::enableIfArithmeticWithReturn<T, T> pow(const T base, const unsigned exponent) noexcept {
  if(exponent == 0) {
    return 1;
  }

  double value = base;

  for(unsigned n = 1; n < exponent; n++) {
    value *= base;
  }

  return value;
}

template<typename T>
PURITY_STRONG constexpr T recPow(const T base, const unsigned exponent) noexcept {
  if(exponent == 0) {
    return 1;
  }

  if(exponent == 1) {
    return base;
  }

  if(exponent % 2 == 0) {
    auto halfProblem = recPow(base, exponent / 2);
    return halfProblem * halfProblem;
  }

  return base * recPow(base, exponent - 1);
}

/*!
 * @brief Integer version that just calls the unsigned power function and inverts the result
 * @warning lots can go wrong here!
 */
template<typename T>
PURITY_STRONG constexpr Traits::enableIfArithmeticWithReturn<T, double> pow(const T base, const int exponent) noexcept {
  if(exponent < 0) {
    return 1.0 / pow(base, static_cast<unsigned>(Temple::Math::abs(exponent)));
  }

  if(exponent == 0) {
    return 1;
  }

  return pow(base, static_cast<unsigned>(exponent));
}

/* Implements Newton's iteration to compute the square root of a positive number
 */
template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> sqrt(const T x) {
  if(x < 0) {
    throw "Square-root domain error: Only real if x >= 0!";
  }

  const T epsilon = std::numeric_limits<T>::epsilon();
  T value = 1;
  T previous = 2;

  while(Temple::Math::abs(previous - value) > epsilon) {
    // store the previous value
    previous = value;

    // compute next iteration
    value = 0.5 * (value + x / value);
  }

  return value;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfIntegralWithReturn<T, T> factorial(const T x) {
  if(x < 0) {
    throw "Factorial domain error!";
  }

  if(x == 0) {
    return 1;
  }

  return x * factorial(x - 1);
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> ln(const T x) {
  unsigned decimalReduction = 0;
  T calcX = x;

  while(abs(calcX) > 10) {
    calcX /= 10;
    decimalReduction += 1;
  }

  // Ensure last division leads to value closer to 1
  if(Math::abs(calcX / 10 - 1) < Math::abs(calcX - 1)) {
    calcX /= 10;
    decimalReduction += 1;
  }

  return Detail::lnSeries(calcX) + decimalReduction * M_LN10;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> log10(const T x) {
  if(x <= 0) {
    throw "Log10 domain error!";
  }

  /* ln(z) = ln(10) * log10(z)
   * -> log10(z) = ln(z) / ln(10)
   */
  return ln(x) / M_LN10;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> log(const T x, const T base) {
  if(x <= 0) {
    throw "Log domain error!";
  }

  /* ln(z) = ln(b) * log_b(z)
   * -> log_b(z) = ln(z) / ln(b)
   */
  return ln(x) / ln(base);
}

/* Implements the infinite series where the derivative is expanded as a binomial
 * series and every term is integrated. Deviates most strongly from std::asin at
 * values very close to the borders. Perhaps it is best to use the approximate
 * form there?
 */
template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> asin(const T x) {
  if(!(-1 <= x && x <= 1)) {
    throw "Inverse sine domain error: only real if -1 < x < 1!";
  }

  if(Temple::Math::abs(x) > 0.90) {
    return (x >= 0) ? Detail::asinApprox(x) : -Detail::asinApprox(-x);
  }

  const T epsilon = std::numeric_limits<T>::epsilon();

  T value = x;
  T upper_factorial = 1;
  T lower_factorial = 1;
  T term = 1;

  for(unsigned n = 1; Temple::Math::abs(term) > epsilon; ++n) {
    upper_factorial *= 2 * (n - 1) + 1;
    lower_factorial *= 2 * n;

    term = (
      upper_factorial / lower_factorial
    ) * pow(x, 2 * n + 1) / (2 * n + 1);

    if(Traits::isnan(term)) {
      break;
    }

    value += term;
  }

  return value;
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> acos(const T x) {
  if(!(-1 <= x && x <= 1)) {
    throw "Inverse cosine domain error: only real if -1 <= x <= 1!";
  }

  return M_PI / 2 - asin(x);
}

template<typename T>
PURITY_STRONG constexpr Traits::enableIfFloatingWithReturn<T, T> atan(const T x) {
  if(!(-M_PI / 2 < x && x < M_PI / 2)) {
    throw "Inverse cosine domain error: only real if -1 < x < 1!";
  }

  return asin(
    x / sqrt(x * x + 1)
  );
}

} // namespace Math
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
