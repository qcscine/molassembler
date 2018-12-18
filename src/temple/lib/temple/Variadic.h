/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Variadic functions
 *
 * Provides functions that can apply to as many Containers as needed.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_VARIADIC_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_VARIADIC_H

#include "temple/constexpr/Optional.h"
#include "temple/constexpr/TupleType.h"
#include "temple/ContainerTraits.h"

#include <algorithm>

namespace temple {

namespace variadic {

namespace detail {

template<typename Container>
std::enable_if_t<
  ::temple::traits::hasSize<Container>::value,
  Optional<std::size_t>
> size(const Container& container) {
  return Optional<std::size_t> {container.size()};
}

template<typename Container>
std::enable_if_t<
  !::temple::traits::hasSize<Container>::value,
  Optional<std::size_t>
> size(const Container& /* container */) {
  return {};
}

} // namespace detail

/**
 * @brief Calculates a lower bound to the size of multiple containers, i.e. sums
 *   up the sizes of all passed containers that have a size member function.
 *
 * @param containers The containers to sum up sizes of
 *
 * @return A lower bound to the amount of elements in all containers combined
 */
template<typename ... Containers>
std::size_t sizeLowerBound(const Containers& ... containers) {
  std::array<Optional<std::size_t>, sizeof...(containers)> sizeOptionals {{
    detail::size(containers)...
  }};

  std::size_t sum = 0;

  for(const auto& sizeOptional : sizeOptionals) {
    if(sizeOptional) {
      sum += sizeOptional.value();
    }
  }

  return sum;
}

namespace detail {

template<typename Vector>
Vector concatenateHelper(Vector& vector) {
  return vector;
}

template<typename Vector, class Container, typename ... Containers>
Vector concatenateHelper(
  Vector& vector,
  const Container& container,
  Containers... containers
) {
  vector.insert(
    std::end(vector),
    std::begin(container),
    std::end(container)
  );

  return concatenateHelper(vector, containers...);
}

} // namespace detail

/*!
 * @brief Concatenate various types of containers together with the same ValueType
 *
 * Requires that each container implements a begin and end iterator.
 */
template<typename... Containers>
auto concatenate(const Containers& ... containers) {
  using ValueTypes = std::tuple<
    traits::getValueType<Containers>...
  >;

  using T = std::tuple_element_t<0, ValueTypes>;

  static_assert(
    TupleType::countType<
      ValueTypes,
      T
    >() == std::tuple_size<ValueTypes>::value,
    "Value types of all containers involved in concatenation must be identical!"
  );

  std::vector<T> concatenated;

  concatenated.reserve(sizeLowerBound(containers...));

  return detail::concatenateHelper(concatenated, containers...);
}

} // namespace variadic

} // namespace temple

#endif
