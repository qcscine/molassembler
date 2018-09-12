#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_VARIADIC_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_VARIADIC_H

#include "temple/constexpr/Optional.h"
#include "temple/constexpr/TupleType.h"
#include "temple/Traits.h"

#include <algorithm>
#include <vector>

/* TODO
 * - variadic size returning an optional<SumType> (would permit variadic
 *   concatenate to preallocate space)
 */

/*!@file
 *
 * @brief Variadic functions
 *
 * Provides functions that can apply to as many Containers as needed.
 */

namespace temple {

namespace variadic {

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

//! Concatenate various types of containers together with the same ValueType
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

  return detail::concatenateHelper(concatenated, containers...);
}

} // namespace variadic

} // namespace temple

#endif
