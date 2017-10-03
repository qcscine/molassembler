#ifndef INCLUDE_CONSTEXPR_MAGIC_MATH_H
#define INCLUDE_CONSTEXPR_MAGIC_MATH_H

#include <cmath>
#include <limits>
#include <type_traits>

/*! @file
 *
 * Provides \c constexpr basic mathematical function implementations, some 
 * logical functions and floating-point comparison helpers.
 */

/* TODO 
 * - Add periodicities of the trigonometric functions
 * - Investigate ill-conditioned quality of inverse trig functions at specific
 *   angles (see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions)
 *
 *   - Patched somewhat with approximation, but still no better than 1e-10 over
 *     the entire domain
 *
 * - Improve power function, it's probably hella awful
 * - Several functions are marked noexcept but throw!
 */
namespace ConstexprMagic {

namespace Math {

namespace traits {

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

} // namespace traits

/* Logic */
//! Template parameter-pack exclusive or of booleans
template<typename ... Bools>
constexpr bool XOR(Bools ... bools);

/* Some very basic math functions for all types */
template<typename T>
inline constexpr T abs(const T& x) noexcept;

template<typename T>
constexpr T max(const T& a, const T& b) noexcept;

template<typename T>
constexpr T min(const T& a, const T& b) noexcept;

/* Floating-point math functions */

// Angle conversions
template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, T> toRadians(const T& inDegrees) noexcept;

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, T> toDegrees(const T& inRadians) noexcept;

/* Comparison helpers, deprecated in favor of identical implementations in
 * FloatingPointComparison.h
 */
template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, bool> isCloseRelative(
  const T& a,
  const T& b,
  const T& relativeTolerance
) __attribute__ ((deprecated));

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, bool> isCloseAbsolute(
  const T& a,
  const T& b,
  const T& absoluteTolerance
) __attribute__ ((deprecated));

// Rounding
template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, int> ceil(const T& value);

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, int> floor(const T& value);

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, int> round(const T& value);

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, T> round(
  const T& value,
  const unsigned& nDigits
);

// Powers
template<typename T>
constexpr T pow(const T& base, const unsigned& exponent) noexcept;

template<typename T>
constexpr T pow(const T& base, const int& exponent) noexcept;

// Sqrt
template<typename T>
constexpr T sqrt(const T& x) noexcept;

// Factorial
template<typename T>
constexpr traits::enableIfIntegralWithReturn<T, T> factorial(const T& x) noexcept;

// Logarithms
template<typename T>
constexpr T ln(const T& x);

template<typename T>
constexpr T log10(const T& x);

template<typename T>
constexpr T log(const T& x, const T& base);

// Inverse trigonometry
template<typename T>
constexpr T asin(const T& x) noexcept;

template<typename T>
constexpr T acos(const T& x);

template<typename T>
constexpr T atan(const T& x);


/* Implementations begin here ------------------------------------------------*/

namespace detail { // Implementation helpers

// Specialization of TPPSum for empty parameter pack
constexpr unsigned TPPSum() {
  return 0;
}

// Template parameter pack sum (needed for XOR function)
template<typename T1, typename... T>
constexpr unsigned TPPSum(T1 a, T ... pack) {
  return a + TPPSum(pack ...);
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, bool> isCloseRelativeOrAbsolute(
  const T& a,
  const T& b,
  const T& relativeTolerance,
  const T& absoluteTolerance
) {
  if(!(
    a != std::numeric_limits<T>::infinity()
    && a != - std::numeric_limits<T>::infinity()
    && b != std::numeric_limits<T>::infinity()
    && b != - std::numeric_limits<T>::infinity()
    && a != std::numeric_limits<T>::quiet_NaN()
    && b != std::numeric_limits<T>::quiet_NaN()
    && a != std::numeric_limits<T>::signaling_NaN()
    && b != std::numeric_limits<T>::signaling_NaN()
  )) {
    throw "isCloseRelativeOrAsbolute cannot handle infinities or NaNs!";
  }
  if(!(relativeTolerance >= 0 && absoluteTolerance >= 0)) {
    throw "isCloseRelativeOrAbsolute: One of either tolerances "
      "needs to be above zero!";
  }

  return(
    abs(abs(a) - abs(b))
    <= max(
      relativeTolerance * max(
        abs(a),
        abs(b)
      ),
      absoluteTolerance
    )
  );
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, int> roundImpl(const T& value) {
  if(value * 10 - floor(value) * 10 >= 5) {
    return ceil(value);
  }

  return floor(value);
}

/* Based on series expansion of ln x:
 * 
 *   y(x) = (x - 1) / (x + 1)
 *   ln x = 2 [ y + y^3/3 + y^5/5 + ...]
 *
 *   for Re x ≥ 0 and x ≠ 0
 */
template<typename T>
constexpr T lnSeries(const T& x) {
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

  for(unsigned n = 3; ConstexprMagic::Math::abs(previous - value) > epsilon; n += 2) {
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
constexpr T asinApprox(const T& x) noexcept {
  if(!(0 < x && x < 1)) {
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

} // namespace detail




template<typename ... Bools>
constexpr bool XOR(Bools ... bools) {
  return detail::TPPSum(bools ...) == 1;
}

template<typename T>
inline constexpr T abs(const T& x) noexcept {
  return (x >= 0) ? x : -x;
}

template<typename T>
constexpr T max(const T& a, const T& b) noexcept {
  return (a > b) ? a : b;
}

template<typename T>
constexpr T min(const T& a, const T& b) noexcept {
  return (a < b) ? a : b;
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, T> toRadians(const T& inDegrees) noexcept {
  return M_PI * inDegrees / 180;
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, T> toDegrees(const T& inRadians) noexcept {
  return 180 * inRadians / M_PI;
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, bool> isCloseRelative(
  const T& a,
  const T& b,
  const T& relativeTolerance
) {
  return detail::isCloseRelativeOrAbsolute(
    a,
    b,
    relativeTolerance,
    T {0}
  );
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, bool> isCloseAbsolute(
  const T& a,
  const T& b,
  const T& absoluteTolerance
) {
  return detail::isCloseRelativeOrAbsolute(
    a,
    b,
    T {0},
    absoluteTolerance
  );
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, int> ceil(const T& value) {
  // Truncate to an int
  const int truncated = static_cast<int>(value);

  if(truncated < value) {
    return truncated + 1;
  }

  return truncated;
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, int> floor(const T& value) {
  // Truncate to an int
  const int truncated = static_cast<int>(value);

  if(truncated > value) {
    return truncated - 1;
  }

  return truncated;
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, int> round(const T& value) {
  if(value < 0) {
    return -detail::roundImpl(-value);
  }

  return detail::roundImpl(value);
}

template<typename T>
constexpr traits::enableIfFloatingWithReturn<T, T> round(
  const T& value,
  const unsigned& nDigits
) {
  if(value == 0) {
    return 0;
  }

  const T d = ceil(log10(abs(value)));
  const int power = nDigits - static_cast<int>(d);
  const T magnitude = pow(10, power);

  // This isn't an integer division since magnitude is a floating-point type
  return round(value * magnitude) / magnitude;
}

// Really weak first implementation
template<typename T>
constexpr T pow(const T& base, const unsigned& exponent) noexcept {
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
constexpr T recPow(const T& base, const unsigned& exponent) noexcept {
  if(exponent == 1) {
    return base;
  }

  if(exponent % 2 == 0) {
    auto halfProblem = recPow(base, exponent / 2);
    return halfProblem * halfProblem;
  }

  return base * recPow(base, exponent - 1);
}

/* Integer version just calls the unsigned power function
 * TODO lots can go wrong here!
 */
template<typename T>
constexpr T pow(const T& base, const int& exponent) noexcept {
  if(exponent < 0) {
    return 1 / pow(base, static_cast<unsigned>(ConstexprMagic::Math::abs(exponent)));
  } 

  if(exponent == 0) {
    return 1;
  }
  
  return pow(base, static_cast<unsigned>(exponent));
}

/* Implements Newton's iteration to compute the square root of a positive number
 */
template<typename T>
constexpr T sqrt(const T& x) noexcept {
  if(x < 0) {
    throw "Square-root domain error: Only real if x >= 0!";
  }

  const T epsilon = std::numeric_limits<T>::epsilon();
  T value = 1;
  T previous = 2;

  while(ConstexprMagic::Math::abs(previous - value) > epsilon) {
    // store the previous value
    previous = value;

    // compute next iteration
    value = 0.5 * (value + x / value);
  }

  return value;
}

template<typename T>
constexpr traits::enableIfIntegralWithReturn<T, T> factorial(const T& x) noexcept {
  if(x <= 0) {
    throw "Factorial domain error!";
  }

  if(x == 1) {
    return 1;
  } else {
    return x * factorial(x - 1);
  }
}

template<typename T>
constexpr T ln(const T& x) {
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

  return detail::lnSeries(calcX) + decimalReduction * M_LN10;
}

template<typename T>
constexpr T log10(const T& x) {
  if(x <= 0) {
    throw "Log10 domain error!";
  }

  /* ln(z) = ln(10) * log10(z)
   * -> log10(z) = ln(z) / ln(10)
   */
  return ln(x) / M_LN10;
}

template<typename T>
constexpr T log(const T& x, const T& base) {
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
constexpr T asin(const T& x) noexcept {
  if(!(-1 < x && x < 1)) {
    throw "Inverse sine domain error: only real if -1 < x < 1!";
  }

  if(ConstexprMagic::Math::abs(x) > 0.92) {
    return (x > 0) ? detail::asinApprox(x) : -detail::asinApprox(-x);
  }

  const T epsilon = std::numeric_limits<T>::epsilon();

  T value = x;
  T upper_factorial = 1;
  T lower_factorial = 1;
  T term = 1;

  for(unsigned n = 1; ConstexprMagic::Math::abs(term) > epsilon; ++n) {
    upper_factorial *= 2 * (n - 1) + 1;
    lower_factorial *= 2 * n;

    term = (
      upper_factorial / lower_factorial
    ) * pow(x, 2 * n + 1) / (2 * n + 1);

    if(std::isnan(term)) {
      break;
    }

    value += term;
  }

  return value;
}

template<typename T>
constexpr T acos(const T& x) {
  if(!(-1 < x && x < 1)) {
    throw "Inverse cosine domain error: only real if -1 < x < 1!";
  }

  return M_PI / 2 - asin(x);
}

template<typename T>
constexpr T atan(const T& x) {
  if(!(-M_PI / 2 < x && x < M_PI / 2)) {
    throw "Inverse cosine domain error: only real if -1 < x < 1!";
  }

  return asin(
    x / sqrt(x * x + 1)
  );
}

} // namespace Math

} // namespace ConstexprMagic

#endif
