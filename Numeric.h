#ifndef INCLUDE_TEMPLATE_MAGIC_NUMERIC_H
#define INCLUDE_TEMPLATE_MAGIC_NUMERIC_H

#include "Traits.h"

#include <numeric>
#include <functional>
#include <cmath>

/*! @file
 * 
 * Exposes a small set of functions for working with numbers in containers.
 */

namespace TemplateMagic {

/*! 
 * Composable sum function. Returns the type the container contains, assuming
 * monadic behavior on operator + (value_type + value_type = value_type).
 * Container must implement begin and end members.
 */
template<class ContainerType>
traits::getValueType<ContainerType>
sum(const ContainerType& container) {
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
std::enable_if_t<
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
double average(const ContainerType& container) {
  assert(&(*container.begin()) != &(*container.end()));

  return static_cast<double>(
    sum(container)
  ) / container.size();
}

template<class ContainerType>
double geometricMean(const ContainerType& container) {
  assert(&(*container.begin()) != &(*container.end()));

  using ValueType = traits::getValueType<ContainerType>;

  return std::pow(
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
 * must have operator+ and be convertible to double.
 */
template<class ContainerType>
double stddev(const ContainerType& container) {
  using ValueType = traits::getValueType<ContainerType>;
  assert(&(*container.begin()) != &(*container.end()));

  const auto averageValue = average(container);

  return std::sqrt(
    sum(
      map(
        container,
        [&averageValue](const ValueType& value) -> ValueType {
          double diff = averageValue - value;
          return diff * diff;
        }
      )
    ) / container.size()
  );
}

//! Variant with known average value
template<class ContainerType>
double stddev(
  const ContainerType& container,
  const double& averageValue
) {
  using ValueType = traits::getValueType<ContainerType>;
  assert(&(*container.begin()) != &(*container.end()));

  return std::sqrt(
    sum(
      map(
        container,
        [&averageValue](const ValueType& value) -> ValueType {
          double diff = averageValue - value;
          return diff * diff;
        }
      )
    ) / container.size()
  );
}

/*!
 * Composable min function. Returns the smallest member of any container.
 *
 * Container must implement begin, end iterators. The iterators must be
 * copy-assignable. The contained type must implement operator <.
 */
template<class ContainerType>
auto min(const ContainerType& container) {
  auto smallestIter = container.begin();
  for(auto it = container.begin(); it != container.end(); ++it) {
    if(*it < *smallestIter) {
      smallestIter = it;
    }
  }

  return *smallestIter;
}

/*!
 * Composable max function. Returns the smallest member of any container.
 *
 * Container must implement begin, end iterators. The iterators must be
 * copy-assignable. The contained type must implement operator <.
 */
template<class ContainerType>
auto max(const ContainerType& container) {
  auto largestIter = container.begin();

  for(auto it = container.begin(); it != container.end(); ++it) {
    if(*largestIter < *it) {
      largestIter = it;
    }
  }

  return *largestIter;
}

} // namespace TemplateMagic

#endif
