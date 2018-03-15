#ifndef INCLUDE_TEMPLATE_MAGIC_FUNCTIONAL_H
#define INCLUDE_TEMPLATE_MAGIC_FUNCTIONAL_H

#include "Invoke.h"
#include "AddToContainer.h"

/*! @file
 * Better functional-style composability improvements than Containers.h.
 * Compatible with range adaptors from Adaptors.h, but very incomplete.
 */

namespace temple {

template<
  typename Container,
  std::enable_if_t<traits::hasSize<Container>::value, int> = 0
> auto size(const Container& container) {
  return container.size();
}

template<
  typename Container,
  std::enable_if_t<!traits::hasSize<Container>::value, int> = 0
> auto size(const Container& container) {
  return std::distance(
    container.begin(),
    container.end()
  );
}

/*!
 * Maps the values in a container using a unary function.
 *
 * Requires:
 * - Container must have template parameters of the form:
 *   Container<ValueType, A<ValueType>, B<ValueType>, ...>,
 *   where zero dependent template parameters (A, B, ...) are also acceptable
 * - Any dependent template parameters (A, B, ..) must be instantiable for the
 *   function return type.
 * - Container implements begin and end forward iterators.
 * - Container implements either insert, emplace, push_back or emplace_back
 * - UnaryFunction must be unary and callable with Container's ValueType
 *
 * Besides custom containers that fulfill the required criteria, this function
 * should be valid for the following STL containers:
 *
 * - vector
 * - deque
 * - list, forward_list
 * - set, multiset, unordered_set, unordered_multiset
 *
 * Notably absent: array, map, multimap, unordered_multimap (some are covered by
 * specializations below!)
 *
 */
template<
  class UnaryFunction,
  typename T,
  template<typename, typename...> class Container,
  template<typename> class ... Dependents
> auto map(
  const Container<T, Dependents<T>...>& container,
  UnaryFunction&& function
) {
  using U = decltype(
    invoke(function, *container.begin())
  );

  Container<U, Dependents<U>...> returnContainer;

  for(const auto& element : container) {
    addToContainer(
      returnContainer,
      invoke(function, element)
    );
  }

  return returnContainer;
}

} // namespace temple

#endif
