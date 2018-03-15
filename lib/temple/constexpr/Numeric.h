#ifndef INCLUDE_TEMPLE_CONSTEXPR_NUMERIC_H
#define INCLUDE_TEMPLE_CONSTEXPR_NUMERIC_H

#include "../Traits.h"
#include "Math.h"

#include <numeric>
#include <functional>
#include <cassert>

/*! @file
 * 
 * Exposes a small set of functions for working with numbers in containers.
 *
 * - Basic summation
 * - Kahan summation
 * - Average
 * - Geometric mean
 * - Standard deviation
 * - Minimum and maximum element
 */

namespace temple {

/*! 
 * Composable sum function. Returns the type the container contains, assuming
 * monadic behavior on operator + (value_type + value_type = value_type).
 * Container must implement begin and end members.
 */
template<class ContainerType>
constexpr traits::getValueType<ContainerType> sum(const ContainerType& container) {
  using ValueType = traits::getValueType<ContainerType>;

  return std::accumulate(
    container.begin(),
    container.end(),
    ValueType {0},
    std::plus<ValueType>()
  );
}

/*! 
 * Composable Kahan summation function. Returns the type the container
 * contains, assuming monadic behavior on operator + (value_type + value_type =
 * value_type).  Container must implement begin and end members.
 */
template<class ContainerType>
constexpr std::enable_if_t<
  std::is_floating_point<
    traits::getValueType<ContainerType>
  >::value,
  traits::getValueType<ContainerType>
> kahanSum(const ContainerType& container) {
  using ValueType = traits::getValueType<ContainerType>;

  ValueType counter {0};
  ValueType error {0};

  for(const auto& value: container) {
    ValueType y = value - error;
    ValueType kt = counter + y;

    error = (kt - counter) - y;
    counter = kt;
  }

  return counter;
}

/*!
 * Composable average function. Returns the average (always as a double) of any
 * basic data type. An assert checks whether the function receives an empty
 * container.
 *
 * Container must implement begin, end and size members, the contained type
 * must have operator+ and be convertible to double.
 */
template<class ContainerType>
constexpr std::enable_if_t<
  std::is_floating_point<traits::getValueType<ContainerType>>::value,
  traits::getValueType<ContainerType>
> average(const ContainerType& container) {
  if(&(*container.begin()) == &(*container.end())) {
    throw "Average called on an empty container!";
  }

  return sum(container) / container.size();
}

template<class ContainerType>
constexpr std::enable_if_t<
  !std::is_floating_point<traits::getValueType<ContainerType>>::value,
  double
> average(const ContainerType& container) {
  if(&(*container.begin()) == &(*container.end())) {
    throw "Average called on an empty container!";
  }

  return static_cast<double>(
    sum(container)
  ) / container.size();
}

template<class ContainerType>
constexpr std::enable_if_t<
  std::is_floating_point<traits::getValueType<ContainerType>>::value,
  traits::getValueType<ContainerType>
> geometricMean(const ContainerType& container) {
  if(&(*container.begin()) == &(*container.end())) {
    throw "geometricMean called on an empty container!";
  }

  using ValueType = traits::getValueType<ContainerType>;

  return Math::pow(
    reduce(
      container,
      ValueType {1.0},
      std::multiplies<ValueType>()
    ),
    1.0 / container.size()
  );
}

/*!
 * Container must implement begin, end and size members, the contained type
 * must have operator+.
 */
template<class ContainerType>
constexpr auto stddev(const ContainerType& container) {
  assert(container.begin() != container.end());

  using ValueType = traits::getValueType<ContainerType>;
  assert(&(*container.begin()) != &(*container.end()));

  const auto averageValue = average(container);

  return Math::sqrt(
    sum(
      map(
        container,
        [&averageValue](const ValueType& value) -> ValueType {
          ValueType diff = averageValue - value;
          return diff * diff;
        }
      )
    ) / container.size()
  );
}

//! Variant with known average value
template<class ContainerType, typename FloatingType>
std::enable_if_t<
  std::is_floating_point<FloatingType>::value,
  FloatingType
> stddev(
  const ContainerType& container,
  const FloatingType& averageValue
) {
  assert(container.begin() != container.end());

  using ValueType = traits::getValueType<ContainerType>;

  static_assert(
    std::is_same<ValueType, FloatingType>::value, 
    "The provided average type must match the container vaue type!"
  );

  assert(&(*container.begin()) != &(*container.end()));

  return Math::sqrt(
    sum(
      map(
        container,
        [&averageValue](const ValueType& value) -> ValueType {
          FloatingType diff = averageValue - value;
          return diff * diff;
        }
      )
    ) / container.size()
  );
}

/*! Composable min function. Returns the smallest value of any container.
 *
 * Container must implement begin, end iterators. The iterators must be
 * copy-assignable. The contained type must implement operator <.
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

/*! Composable max function. Returns the smallest value of any container.
 *
 * Container must implement begin, end iterators. The iterators must be
 * copy-assignable. The contained type must implement operator <.
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

} // namespace temple

#endif
