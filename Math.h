#ifndef INCLUDE_CONSTEXPR_MAGIC_MATH_H
#define INCLUDE_CONSTEXPR_MAGIC_MATH_H

#include <cmath>
#include <limits>
#include <type_traits>
#include <cassert>

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

/* Some very basic math functions */
template<typename T>
inline constexpr T abs(const T& x) noexcept {
  return (x >= 0) ? x : -x;
}

namespace detail {

template<typename T>
constexpr T max(const T& a, const T& b) noexcept {
  return (a > b) ? a : b;
}

template<typename T>
constexpr T min(const T& a, const T& b) noexcept {
  return (a < b) ? a : b;
}

template<typename T, typename U>
using lesserOfType = std::conditional_t<
  sizeof(T) <= sizeof(U),
  T,
  U
>;

template<typename T, typename U>
constexpr 
std::enable_if_t<
  std::is_floating_point<T>::value && std::is_floating_point<U>::value,
  bool
> isClose(T a, U b) noexcept {
  return (
    abs(abs(a) - abs(b)) <= (
      std::numeric_limits<lesserOfType<T, U>>::epsilon() 
      * max(abs(a), abs(b))
    )
  );
}

template<typename T, typename U>
constexpr 
std::enable_if_t<
  !std::is_floating_point<T>::value || !std::is_floating_point<U>::value,
  bool
> isClose(T a, U b) noexcept {
  return a == b;
}

template<typename T>
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelativeOrAbsolute(
  const T& a,
  const T& b,
  const T& relativeTolerance,
  const T& absoluteTolerance
) {
  assert(
    a != std::numeric_limits<T>::infinity()
    && a != - std::numeric_limits<T>::infinity()
    && b != std::numeric_limits<T>::infinity()
    && b != - std::numeric_limits<T>::infinity()
    && a != std::numeric_limits<T>::quiet_NaN()
    && b != std::numeric_limits<T>::quiet_NaN()
    && a != std::numeric_limits<T>::signaling_NaN()
    && b != std::numeric_limits<T>::signaling_NaN()
  );
  assert(relativeTolerance >= 0 && absoluteTolerance >= 0);

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
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isApprox(const T& a, const T& b, const T& epsilon) noexcept {
  if(a == b) return true;

  auto absDiff = ConstexprMagic::Math::abs(
    ConstexprMagic::Math::abs(a) 
    - ConstexprMagic::Math::abs(b)
  );
  
  if(a == 0 || b == 0 || absDiff < std::numeric_limits<T>::epsilon()) {
    return absDiff < epsilon * std::numeric_limits<T>::epsilon();
  }
  
  return absDiff / min(
    ConstexprMagic::Math::abs(a) + ConstexprMagic::Math::abs(b), 
    std::numeric_limits<T>::epsilon()
  ) < epsilon;
}

template<typename T>
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  int
> roundImpl(const T& value) {
  if(value * 10 - floor(value) * 10 >= 5) {
    return ceil(value);
  }

  return floor(value);
}

} // namespace detail

template<typename T>
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelative(
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
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseAbsolute(
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
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  int
> ceil(const T& value) {
  // Truncate to an int
  const int truncated = static_cast<int>(value);

  if(truncated < value) {
    return truncated + 1;
  }

  return truncated;
}

template<typename T>
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  int
> floor(const T& value) {
  // Truncate to an int
  const int truncated = static_cast<int>(value);

  if(truncated > value) {
    return truncated - 1;
  }

  return truncated;
}

template<typename T>
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  int
> round(const T& value) {
  if(value < 0) {
    return -detail::roundImpl(-value);
  }

  return detail::roundImpl(value);
}

template<typename T>
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  T
> round(
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

constexpr double toRadians(const double& inDegrees) noexcept {
  return M_PI * inDegrees / 180;
}

constexpr double toDegrees(const double& inRadians) noexcept {
  return 180 * inRadians / M_PI;
}

// Really weak first implementation
constexpr double pow(const double& base, const unsigned& exponent) noexcept {
  double value = base;

  for(unsigned n = 1; n < exponent; n++) {
    value *= base;
  }

  return value;
}

/* Integer version just calls the unsigned power function
 * TODO lots can go wrong here!
 */
constexpr double pow(const double& base, const int& exponent) noexcept {
  if(exponent < 0) {
    return 1 / pow(base, static_cast<unsigned>(ConstexprMagic::Math::abs(exponent)));
  } 

  if(exponent == 0) {
    return 1;
  }
  
  return pow(base, static_cast<unsigned>(exponent));
}


/*template<typename T>
inline constexpr auto floor(const T& x) noexcept -> decltype(std::floor(x)) {
  return (
    (int(x) == x) ? int(x) : (
      (x >= 0.0) ? int(x) : int(x) - 1
    )
  );
}*/

/* Implements Newton's iteration to compute the square root of a positive number
 */
template<typename T>
constexpr T sqrt(const T& x) noexcept {
  assert(
    x >= 0 
    && "Square-root domain error: Only real if x >= 0!"
  );

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

/* Based on series expansion of ln x:
 * 
 *   y(x) = (x - 1) / (x + 1)
 *   ln x = 2 [ y + y^3/3 + y^5/5 + ...]
 *
 *   for Re x ≥ 0 and x ≠ 0
 */
template<typename T>
constexpr T lnSeries(const T& x) {
  assert(x > 0);

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

template<typename T>
constexpr T ln(const T& x) {
  unsigned decimalReduction = 0;
  T calcX = x;

  while(abs(calcX) > 10) {
    calcX /= 10;
    decimalReduction += 1;
  }

  // Ensure last division leads to value closer to 1
  if(std::fabs(calcX / 10 - 1) < std::fabs(calcX - 1)) {
    calcX /= 10;
    decimalReduction += 1;
  }

  return lnSeries(calcX) + decimalReduction * M_LN10;
}

template<typename T>
constexpr T log10(const T& x) {
  assert(x > 0);

  /* ln(z) = ln(10) * log10(z)
   * -> log10(z) = ln(z) / ln(10)
   */
  return ln(x) / M_LN10;
}

/* Implements an approximation to asin over 0 < x < 1 with accuracy > 2e-8
 * from Abramowitz, M., Stegun, I.: Handbook of Mathematical Functions, p. 64,
 * 1964, from http://people.math.sfu.ca/~cbm/aands/abramowitz_and_stegun.pdf
 */
template<typename T>
constexpr T asinApprox(const T& x) noexcept {
  assert(
    0 < x && x < 1 
    && "Asin approximation domain error: only applicable for 0 < x < 1!"
  );

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

/* Implements the infinite series where the derivative is expanded as a binomial
 * series and every term is integrated. Deviates most strongly from std::asin at
 * values very close to the borders. Perhaps it is best to use the approximate
 * form there?
 */
template<typename T>
constexpr T asin(const T& x) noexcept {
  assert(
    -1 < x && x < 1
    && "Inverse sine domain error: only real if -1 < x < 1!"
  );

  if(ConstexprMagic::Math::abs(x) > 0.92) {
    return (x > 0) ? asinApprox(x) : -asinApprox(-x);
  }

  const T epsilon = std::numeric_limits<T>::epsilon();

  T value = x;
  T upper_factorial = 1;
  T lower_factorial = 1;
  T term = 1;

  for(unsigned n = 1; ConstexprMagic::Math::abs(term) > epsilon; n++) {
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
constexpr T acos(const T& x) noexcept {
  if(!(-1 < x && x < 1)) {
    throw "Inverse cosine domain error: only real if -1 < x < 1!";
  }

  return M_PI / 2 - asin(x);
}

} // namespace Math

} // namespace ConstexprMagic

#endif
