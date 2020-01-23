/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Helpers for comparing floating point values
 *
 * Some floating point comparison helpers for direct equality testing on the
 * basis of either relative or absolute tolerance. Two comparator classes for
 * better legibility of common comparison operators in expanded equality
 * contexts.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLTE_CONSTEXPR_FP_COMPARISON_H
#define INCLUDE_MOLASSEMBLER_TEMPLTE_CONSTEXPR_FP_COMPARISON_H

#include "temple/constexpr/Math.h"

#include <tuple>
#include <cassert>

namespace Scine {
namespace temple {

//! @brief Floating-point comparison helpers
namespace floating {

// Floating point comparison helpers
template<typename T>
constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelative(
  T a,
  T b,
  T relativeTolerance
);

template<typename T>
constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseAbsolute(
  T a,
  T b,
  T absoluteTolerance
);


// Implementation details
namespace detail {

template<typename T>
constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelativeOrAbsolute(
  T a,
  T b,
  T relativeTolerance,
  T absoluteTolerance
);

template<typename T>
PURITY_STRONG constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelativeOrAbsolute(
  const T a,
  const T b,
  const T relativeTolerance,
  const T absoluteTolerance
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
    throw "isCloseRelativeOrAbsolute cannot handle infinities or NaNs!";
  }

  if(!(relativeTolerance >= 0 && absoluteTolerance >= 0)) {
    throw "isCloseRelativeOrAbsolute: One of either tolerances "
      "needs to be above zero!";
  }

  return(
    Math::abs(a - b)
    <= Math::max(
      relativeTolerance * Math::max(
        Math::abs(a),
        Math::abs(b)
      ),
      absoluteTolerance
    )
  );
}

} // namespace detail

template<typename T>
PURITY_STRONG constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelative(
  const T a,
  const T b,
  const T relativeTolerance
) {
  return detail::isCloseRelativeOrAbsolute(
    a,
    b,
    relativeTolerance,
    T {0}
  );
}

template<typename T>
PURITY_STRONG constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseAbsolute(
  const T a,
  const T b,
  const T absoluteTolerance
) {
  return detail::isCloseRelativeOrAbsolute(
    a,
    b,
    T {0},
    absoluteTolerance
  );
}

// Comparators exploiting relative and absolute expanded definitions of equality

template<typename T>
class ExpandedAbsoluteEqualityComparator {
private:
  const T _absoluteTolerance;

public:
  constexpr ExpandedAbsoluteEqualityComparator(const T absoluteTolerance)
    : _absoluteTolerance(absoluteTolerance)
  {
    assert(
      absoluteTolerance > 0
      && absoluteTolerance != std::numeric_limits<T>::infinity()
      && absoluteTolerance != std::numeric_limits<T>::quiet_NaN()
      && absoluteTolerance != std::numeric_limits<T>::signaling_NaN()
    );
  }

  PURITY_STRONG constexpr bool isLessThan(const T a, const T b) const noexcept {
    return a < (b - _absoluteTolerance);
  }

  PURITY_STRONG constexpr bool isMoreThan(const T a, const T b) const noexcept {
    return a > (b + _absoluteTolerance);
  }

  PURITY_STRONG constexpr bool isLessOrEqual(const T a, const T b) const noexcept {
    return a < (b + _absoluteTolerance);
  }

  PURITY_STRONG constexpr bool isMoreOrEqual(const T a, const T b) const noexcept {
    return a > (b - _absoluteTolerance);
  }

  PURITY_STRONG constexpr bool isEqual(const T a, const T b) const noexcept {
    return Math::abs(a - b) <= _absoluteTolerance;
  }

  PURITY_STRONG constexpr bool isUnequal(const T a, const T b) const noexcept {
    return !isEqual(a, b);
  }

  //! Function call operator compares equality
  PURITY_STRONG constexpr bool operator () (const T a, const T b) const noexcept {
    return isEqual(a, b);
  }
};

template<typename T>
class ExpandedRelativeEqualityComparator {
private:
  const T _relativeTolerance;

public:
  constexpr ExpandedRelativeEqualityComparator(const T relativeTolerance)
    : _relativeTolerance(relativeTolerance)
  {
    assert(relativeTolerance > 0);
  }

  PURITY_STRONG constexpr bool isLessThan(const T a, const T b) const {
    return (
      (a < b) && !detail::isCloseRelativeOrAbsolute(
        a,
        b,
        _relativeTolerance,
        T {0}
      )
    );
  }

  PURITY_STRONG constexpr bool isMoreThan(const T a, const T b) const {
    return (
      (a > b) && !detail::isCloseRelativeOrAbsolute(
        a,
        b,
        _relativeTolerance,
        T {0}
      )
    );
  }

  PURITY_STRONG constexpr bool isLessOrEqual(const T a, const T b) const {
    return (
      (a < b) || detail::isCloseRelativeOrAbsolute(
        a,
        b,
        _relativeTolerance,
        T {0}
      )
    );
  }

  PURITY_STRONG constexpr bool isMoreOrEqual(const T a, const T b) const {
    return (
      (a > b) || detail::isCloseRelativeOrAbsolute(
        a,
        b,
        _relativeTolerance,
        T {0}
      )
    );
  }

  PURITY_STRONG constexpr bool isEqual(const T a, const T b) const {
    return detail::isCloseRelativeOrAbsolute(
      a,
      b,
      _relativeTolerance,
      T {0}
    );
  }

  PURITY_STRONG constexpr bool isUnequal(const T a, const T b) const {
    return !detail::isCloseRelativeOrAbsolute(
      a,
      b,
      _relativeTolerance,
      T {0}
    );
  }

  //! Function call operator compares equality
  PURITY_STRONG constexpr bool operator () (const T a, const T b) const noexcept {
    return isEqual(a, b);
  }
};

} // namespace floating
} // namespace temple
} // namespace Scine

#endif
