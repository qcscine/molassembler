/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Stubs to work with numeric data
 *
 * Exposes a small set of functions for working with numbers in containers.
 *
 * - Basic summation
 * - Average
 * - Geometric mean
 * - Standard deviation
 * - Minimum and maximum element
 */

#ifndef INCLUDE_TEMPLE_CONSTEXPR_NUMERIC_H
#define INCLUDE_TEMPLE_CONSTEXPR_NUMERIC_H

#include "Molassembler/Temple/Traits.h"
#include "Molassembler/Temple/constexpr/Math.h"

#include <numeric>
#include <functional>
#include <cassert>

namespace Scine {
namespace Molassembler {
namespace Temple {

/*! @brief Summation with zero-initialization
 *
 * Composable sum function. Returns the type the container contains, assuming
 * monadic behavior on operator + (value_type + value_type = value_type).
 * Container must implement begin and end members.
 *
 * @complexity{@math{\Theta(N)}}
 */
template<class ContainerType>
constexpr Traits::getValueType<ContainerType> sum(const ContainerType& container) {
  using ValueType = Traits::getValueType<ContainerType>;

  return std::accumulate(
    std::begin(container),
    std::end(container),
    ValueType {0},
    std::plus<ValueType>()
  );
}

/*!
 * Composable average function. Returns the average (always as a double) of any
 * basic data type. An assert checks whether the function receives an empty
 * container.
 *
 * Container must implement begin, end and size members, the contained type
 * must have operator+ and be convertible to double.
 *
 * @complexity{@math{\Theta(N)}}
 */
template<class ContainerType>
constexpr std::enable_if_t<
  std::is_floating_point<Traits::getValueType<ContainerType>>::value,
  Traits::getValueType<ContainerType>
> average(const ContainerType& container) {
  if(container.begin() == container.end()) {
    throw "Average called on an empty container!";
  }

  return sum(container) / container.size();
}

//! @overload
template<class ContainerType>
constexpr std::enable_if_t<
  !std::is_floating_point<Traits::getValueType<ContainerType>>::value,
  double
> average(const ContainerType& container) {
  if(container.begin() == container.end()) {
    throw "Average called on an empty container!";
  }

  return static_cast<double>(
    sum(container)
  ) / container.size();
}

/**
 * @brief Geometric average of all values in a container
 *
 * @param container
 *
 * @complexity{@math{\Theta(N)}}
 */
template<class ContainerType>
constexpr std::enable_if_t<
  std::is_floating_point<Traits::getValueType<ContainerType>>::value,
  Traits::getValueType<ContainerType>
> geometricMean(const ContainerType& container) {
  if(container.begin() == container.end()) {
    throw "geometricMean called on an empty container!";
  }

  using ValueType = Traits::getValueType<ContainerType>;

  return Math::pow(
    reduce(
      container,
      ValueType {1.0},
      std::multiplies<ValueType>()
    ),
    1.0 / container.size()
  );
}

/*! @brief Calculate the standard deviation of a container with a known average
 *
 * @complexity{@math{\Theta(N)}}
 */
template<class ContainerType, typename FloatingType>
constexpr std::enable_if_t<
  std::is_floating_point<FloatingType>::value,
  FloatingType
> stddev(
  const ContainerType& container,
  const FloatingType averageValue
) {
  assert(container.begin() != container.end());

  FloatingType sum = 0;

  for(const auto value : container) {
    auto diff = static_cast<FloatingType>(value) - averageValue;
    sum += diff * diff;
  }

  return Math::sqrt(sum / container.size());
}


/*! Calculates standard deviation without an existing average value
 *
 * Container must implement begin, end and size members, the contained type
 * must have operator+.
 *
 * @complexity{@math{\Theta(2N)}}
 */
template<class ContainerType>
constexpr auto stddev(const ContainerType& container) {
  assert(container.begin() != container.end());

  return stddev(container, average(container));
}

/*! @brief Composable min_element function. Returns the smallest value of any
 *   container.
 *
 * Container must implement begin, end iterators. The iterators must be
 * copy-assignable. The contained type must implement operator <.
 *
 * @complexity{@math{\Theta(N)}}
 * @todo rename min_element
 */
template<class ContainerType>
constexpr auto min(const ContainerType& container) {
  if(container.begin() == container.end()) {
    throw "Min called on empty container!";
  }

  auto smallestIter = container.begin();
  for(auto it = container.begin(); it != container.end(); ++it) {
    if(*it < *smallestIter) {
      smallestIter = it;
    }
  }

  return *smallestIter;
}

/*! @brief Composable max function. Returns the smallest value of any container.
 *
 * Container must implement begin, end iterators. The iterators must be
 * copy-assignable. The contained type must implement operator <.
 *
 * @complexity{@math{\Theta(N)}}
 * @todo rename max_element
 */
template<class ContainerType>
constexpr auto max(const ContainerType& container) {
  if(container.begin() == container.end()) {
    throw "Max called on empty container!";
  }

  auto largestIter = container.begin();

  for(auto it = container.begin(); it != container.end(); ++it) {
    if(*largestIter < *it) {
      largestIter = it;
    }
  }

  return *largestIter;
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
