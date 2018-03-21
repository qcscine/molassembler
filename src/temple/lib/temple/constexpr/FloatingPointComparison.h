#ifndef INCLUDE_CONSTEXPR_MAGIC_FLOATING_POINT_COMPARISON_H
#define INCLUDE_CONSTEXPR_MAGIC_FLOATING_POINT_COMPARISON_H

#include "Math.h"

#include <tuple>
#include <cassert>

/*! @file
 *
 * Some floating point comparison helpers for direct equality testing on the
 * basis of either relative or absolute tolerance. Two comparator classes for
 * better legibility of common comparison operators in expanded equality
 * contexts.
 */

namespace temple {

namespace floating {

// Floating point comparison helpers
template<typename T>
constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelative(
  const T& a,
  const T& b,
  const T& relativeTolerance
);

template<typename T>
constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseAbsolute(
  const T& a,
  const T& b,
  const T& absoluteTolerance
);


// Implementation details
namespace detail {

template<typename T>
constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelativeOrAbsolute(
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
    throw "isCloseRelativeOrAbsolute cannot handle infinities or NaNs!";
  }

  if(!(relativeTolerance >= 0 && absoluteTolerance >= 0)) {
    throw "isCloseRelativeOrAbsolute: One of either tolerances "
      "needs to be above zero!";
  }

  return(
    Math::abs(Math::abs(a) - Math::abs(b))
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
constexpr std::enable_if_t<
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
constexpr std::enable_if_t<
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

// Comparators exploiting relative and absolute expanded definitions of equality

template<typename T>
class ExpandedAbsoluteEqualityComparator {
private:
  const T _absoluteTolerance;

public:
  constexpr ExpandedAbsoluteEqualityComparator(const T& absoluteTolerance)
    : _absoluteTolerance(absoluteTolerance)
  {
    assert(
      absoluteTolerance > 0
      && absoluteTolerance != std::numeric_limits<T>::infinity()
      && absoluteTolerance != std::numeric_limits<T>::quiet_NaN()
      && absoluteTolerance != std::numeric_limits<T>::signaling_NaN()
    );
  }

  constexpr bool isLessThan(const T& a, const T& b) const noexcept {
    return a < (b - _absoluteTolerance);
  }

  constexpr bool isMoreThan(const T& a, const T& b) const noexcept {
    return a > (b + _absoluteTolerance);
  }

  constexpr bool isLessOrEqual(const T& a, const T& b) const noexcept {
    return a < (b + _absoluteTolerance);
  }

  constexpr bool isMoreOrEqual(const T& a, const T& b) const noexcept {
    return a > (b - _absoluteTolerance);
  }

  constexpr bool isEqual(const T& a, const T& b) const noexcept {
    return Math::abs(Math::abs(a) - Math::abs(b)) <= _absoluteTolerance;
  }

  constexpr bool isUnequal(const T& a, const T& b) const noexcept {
    return !isEqual(a, b);
  }

  //! Function call operator compares equality
  constexpr bool operator () (const T& a, const T& b) const noexcept {
    return isEqual(a, b);
  }
};

template<typename T>
class ExpandedRelativeEqualityComparator {
private:
  const T _relativeTolerance;

public:
  constexpr ExpandedRelativeEqualityComparator(const T& relativeTolerance)
    : _relativeTolerance(relativeTolerance)
  {
    assert(relativeTolerance > 0);
  }

  constexpr bool isLessThan(const T& a, const T& b) const {
    return (
      (a < b) && !detail::isCloseRelativeOrAbsolute(
        a,
        b,
        _relativeTolerance,
        T {0}
      ) 
    );
  }

  constexpr bool isMoreThan(const T& a, const T& b) const {
    return (
      (a > b) && !detail::isCloseRelativeOrAbsolute(
        a,
        b,
        _relativeTolerance,
        T {0}
      ) 
    );
  }

  constexpr bool isLessOrEqual(const T& a, const T& b) const {
    return (
      (a < b) || detail::isCloseRelativeOrAbsolute(
        a,
        b,
        _relativeTolerance,
        T {0}
      )
    );
  }

  constexpr bool isMoreOrEqual(const T& a, const T& b) const {
    return (
      (a > b) || detail::isCloseRelativeOrAbsolute(
        a,
        b,
        _relativeTolerance,
        T {0}
      )
    );
  }

  constexpr bool isEqual(const T& a, const T& b) const {
    return detail::isCloseRelativeOrAbsolute(
      a,
      b,
      _relativeTolerance,
      T {0}
    );
  }

  constexpr bool isUnequal(const T& a, const T& b) const {
    return !detail::isCloseRelativeOrAbsolute(
      a,
      b,
      _relativeTolerance,
      T {0}
    );
  }

  //! Function call operator compares equality
  constexpr bool operator () (const T& a, const T& b) const noexcept {
    return isEqual(a, b);
  }
};

} // namespace floating

} // namespace temple

#endif
