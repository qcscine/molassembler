/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Helpers for comparing floating point values
 *
 * Some floating point comparison helpers for direct equality testing on the
 * basis of either relative or absolute tolerance. Two comparator classes for
 * better legibility of common comparison operators in expanded equality
 * contexts.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLTE_CONSTEXPR_FP_COMPARISON_H
#define INCLUDE_MOLASSEMBLER_TEMPLTE_CONSTEXPR_FP_COMPARISON_H

#include "Molassembler/Temple/constexpr/Math.h"

#include <tuple>
#include <cassert>

namespace Scine {
namespace Molassembler {
namespace Temple {

//! @brief Floating-point comparison helpers
namespace Floating {

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
namespace Detail {

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
    && a == a // not NaN
    && b == b // not NaN
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

} // namespace Detail

template<typename T>
PURITY_STRONG constexpr std::enable_if_t<
  std::is_floating_point<T>::value,
  bool
> isCloseRelative(
  const T a,
  const T b,
  const T relativeTolerance
) {
  return Detail::isCloseRelativeOrAbsolute(
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
  return Detail::isCloseRelativeOrAbsolute(
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
  const T absoluteTolerance_;

public:
  constexpr ExpandedAbsoluteEqualityComparator(const T absoluteTolerance)
    : absoluteTolerance_(absoluteTolerance)
  {
    assert(
      absoluteTolerance > 0
      && absoluteTolerance != std::numeric_limits<T>::infinity()
      && absoluteTolerance == absoluteTolerance // not NaN
    );
  }

  constexpr bool isLessThan(const T a, const T b) const noexcept {
    return a < (b - absoluteTolerance_);
  }

  constexpr bool isMoreThan(const T a, const T b) const noexcept {
    return a > (b + absoluteTolerance_);
  }

  constexpr bool isLessOrEqual(const T a, const T b) const noexcept {
    return a < (b + absoluteTolerance_);
  }

  constexpr bool isMoreOrEqual(const T a, const T b) const noexcept {
    return a > (b - absoluteTolerance_);
  }

  constexpr bool isEqual(const T a, const T b) const noexcept {
    return Math::abs(a - b) <= absoluteTolerance_;
  }

  constexpr bool isUnequal(const T a, const T b) const noexcept {
    return !isEqual(a, b);
  }

  //! Function call operator compares equality
  constexpr bool operator () (const T a, const T b) const noexcept {
    return isEqual(a, b);
  }
};

template<typename T>
class ExpandedRelativeEqualityComparator {
private:
  const T relativeTolerance_;

public:
  constexpr ExpandedRelativeEqualityComparator(const T relativeTolerance)
    : relativeTolerance_(relativeTolerance)
  {
    assert(relativeTolerance > 0);
  }

  constexpr bool isLessThan(const T a, const T b) const {
    return (
      (a < b) && !Detail::isCloseRelativeOrAbsolute(
        a,
        b,
        relativeTolerance_,
        T {0}
      )
    );
  }

  constexpr bool isMoreThan(const T a, const T b) const {
    return (
      (a > b) && !Detail::isCloseRelativeOrAbsolute(
        a,
        b,
        relativeTolerance_,
        T {0}
      )
    );
  }

  constexpr bool isLessOrEqual(const T a, const T b) const {
    return (
      (a < b) || Detail::isCloseRelativeOrAbsolute(
        a,
        b,
        relativeTolerance_,
        T {0}
      )
    );
  }

  constexpr bool isMoreOrEqual(const T a, const T b) const {
    return (
      (a > b) || Detail::isCloseRelativeOrAbsolute(
        a,
        b,
        relativeTolerance_,
        T {0}
      )
    );
  }

  constexpr bool isEqual(const T a, const T b) const {
    return Detail::isCloseRelativeOrAbsolute(
      a,
      b,
      relativeTolerance_,
      T {0}
    );
  }

  constexpr bool isUnequal(const T a, const T b) const {
    return !Detail::isCloseRelativeOrAbsolute(
      a,
      b,
      relativeTolerance_,
      T {0}
    );
  }

  //! Function call operator compares equality
  constexpr bool operator () (const T a, const T b) const noexcept {
    return isEqual(a, b);
  }
};

} // namespace Floating
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
