#ifndef INCLUDE_CONSTEXPR_MAGIC_MATH_H
#define INCLUDE_CONSTEXPR_MAGIC_MATH_H

#include <cmath>
#include <limits>
#include <type_traits>

/* TODO 
 * - Add periodicities of the trigonometric functions
 * - Investigate ill-conditioned quality of inverse trig functions at specific
 *   angles (see https://en.wikipedia.org/wiki/Inverse_trigonometric_functions)
 *
 *   - Patched somewhat with approximation, but still no better than 1e-10 over
 *     the entire domain
 *
 * - Improve power function, it's probably hella awful
 */
namespace ConstexprMagic {

namespace Math {

// Not a whole lot you can do wrong here
template<typename T>
inline constexpr T abs(const T& x) {
  return (x >= 0) ? x : -x;
}

namespace detail {

template<typename T>
constexpr bool max(const T& a, const T& b) {
  return (a > b) ? a : b;
}

template<typename T>
constexpr bool min(const T& a, const T& b) {
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
> isClose(T a, U b) {
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
> isClose(T a, U b) {
  return a == b;
}

template<typename T>
constexpr
std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isApprox(const T& a, const T& b, const T& epsilon) {
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

} // namespace detail

// Really weak first implementation
constexpr double pow(const double& base, const unsigned& exponent) {
  double value = base;

  for(unsigned n = 1; n < exponent; n++) {
    value *= base;
  }

  return value;
}

/* Integer version just calls the unsigned power function
 * TODO lots can go wrong here!
 */
constexpr double pow(const double& base, const int& exponent) {
  if(exponent < 0) {
    return 1 / pow(base, static_cast<unsigned>(ConstexprMagic::Math::abs(exponent)));
  }

  else return pow(base, static_cast<unsigned>(exponent));
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
template<typename T>
constexpr T sqrt(const T& x) {
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

/* Implements an approximation to asin over 0 < x < 1 with accuracy > 2e-8
 * from Abramowitz, M., Stegun, I.: Handbook of Mathematical Functions, p. 64,
 * 1964, from http://people.math.sfu.ca/~cbm/aands/abramowitz_and_stegun.pdf
 */
template<typename T>
constexpr T asinApprox(const T& x) {
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

/* Implements the infinite series where the derivative is expanded as a binomial
 * series and every term is integrated. Deviates most strongly from std::asin at
 * values very close to the borders. Perhaps it is best to use the approximate
 * form there?
 */
template<typename T>
constexpr T asin(const T& x) {
  if(!(-1 < x && x < 1)) {
    throw "Inverse sine domain error: only real if -1 < x < 1!";
  }

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
constexpr T acos(const T& x) {
  if(!(-1 < x && x < 1)) {
    throw "Inverse cosine domain error: only real if -1 < x < 1!";
  }

  return M_PI / 2 - asin(x);
}

} // namespace Math

} // namespace ConstexprMagic

#endif
