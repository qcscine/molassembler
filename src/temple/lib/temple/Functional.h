// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_FUNCTIONAL_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_FUNCTIONAL_H

#include "temple/Invoke.h"
#include "temple/AddToContainer.h"

#include <algorithm>
#include <numeric>

/*! @file
 *
 * @brief Functional-style container-related algorithms
 *
 * Supplies a large number of functional-style algorithm definitions. Most are
 * STL algorithm shortcuts but many are composable with the range adaptors
 * provided in Adaptors/ using the function call template functions in
 * Invoke.h.
 */

namespace temple {

template<
  class Container,
  std::enable_if_t<traits::hasSize<Container>::value, int> = 0
> auto size(const Container& container) {
  return container.size();
}

template<
  class Container,
  std::enable_if_t<!traits::hasSize<Container>::value, int> = 0
> auto size(const Container& container) {
  return std::distance(
    std::begin(container),
    std::end(container)
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
> auto map_stl(
  const Container<T, Dependents<T>...>& container,
  UnaryFunction&& function
) {
  using U = decltype(
    invoke(function, *std::begin(container))
  );

  Container<U, Dependents<U>...> returnContainer;
  reserveIfPossible(returnContainer, container);

  for(const auto& element : container) {
    addToContainer(
      returnContainer,
      invoke(function, element)
    );
  }

  return returnContainer;
}

template<
  template<typename, std::size_t> class ArrayType,
  typename T,
  std::size_t size,
  class UnaryFunction
> auto map_stl(
  const ArrayType<T, size>& container,
  UnaryFunction&& function
) {
  using U = decltype(
    invoke(function, *std::begin(container))
  );

  ArrayType<U, size> returnContainer;

  for(unsigned i = 0; i < size; ++i) {
    returnContainer[i] = invoke(function, container[i]);
  }

  return returnContainer;
}

template<class Container, class UnaryFunction>
auto map(
  const Container& container,
  UnaryFunction&& function
) {
  using U = decltype(
    invoke(function, *std::begin(container))
  );

  std::vector<U> returnContainer;
  reserveIfPossible(returnContainer, container);

  for(const auto& value : container) {
    returnContainer.emplace_back(
      invoke(function, value)
    );
  }

  return returnContainer;
}

template<class Container, class Callable>
void forEach(
  const Container& container,
  Callable&& callable
) {
  for(const auto& value : container) {
    invoke(callable, value);
  }
}

template<
  class Container,
  typename ComparisonFunction,
  typename MappingFunction
>
typename Container::const_iterator select(
  const Container& container,
  ComparisonFunction&& comparator,
  MappingFunction&& mapFunction
) {
  auto selection = std::begin(container);
  auto selectionValue = mapFunction(*selection);

  auto iter = std::begin(container);

  while(iter != std::end(container)) {
    auto currentValue = mapFunction(*iter);

    // Replace if 'better' according to the comparator
    if(comparator(currentValue, selectionValue)) {
      selection = iter;
      selectionValue = currentValue;
    }

    ++iter;
  }

  return selection;
}

template<
  class Container,
  typename T,
  class BinaryFunction
> T accumulate(
  const Container& container,
  T init,
  BinaryFunction&& reductionFunction
) {

  for(const auto& value: container) {
    // Call variadic invoke (no tuple unpacking!)
    init = invoke(reductionFunction, init, value);
  }

  return init;
}

template<class Container, class UnaryPredicate>
bool all_of(const Container& container, UnaryPredicate&& predicate) {
  for(const auto& element : container) {
    if(!invoke(predicate, element)) {
      return false;
    }
  }

  return true;
}

template<class Container, class UnaryPredicate>
bool any_of(const Container& container, UnaryPredicate&& predicate) {
  for(const auto& element : container) {
    if(invoke(predicate, element)) {
      return true;
    }
  }

  return false;
}

namespace inplace {

//! Call std::sort on a container
template<class Container>
void sort(Container& container) {
  std::sort(
    std::begin(container),
    std::end(container)
  );
}

//! Call std::sort with a custom comparator on a container
template<class Container, typename Comparator>
void sort(Container& container, Comparator&& comparator) {
  std::sort(
    std::begin(container),
    std::end(container),
    std::forward<Comparator>(comparator)
  );
}

template<class Container, class UnaryFunction>
void transform(Container& container, UnaryFunction&& function) {
  for(auto& value : container) {
    value = function(value);
  }
}

//! Call std::reverse
template<class Container>
void reverse(Container& container) {
  std::reverse(
    std::begin(container),
    std::end(container)
  );
}

//! Apply erase-remove idiom to container
template<class Container, typename T>
void remove(
  Container& container,
  const T& value
) {
  container.erase(
    std::remove(
      std::begin(container),
      std::end(container),
      value
    ),
    std::end(container)
  );
}

//! Apply erase-remove idiom to container with a predicate function
template<class Container, class UnaryFunction>
void remove_if(
  Container& container,
  UnaryFunction&& predicate
) {
  // TODO if container's value type is a tuple and predicate cannot be called
  // with it directly, we have to generate a callable object that forwards an
  // unpacked call!

  container.erase(
    std::remove_if(
      std::begin(container),
      std::end(container),
      std::forward<UnaryFunction>(predicate)
    ),
    std::end(container)
  );
}

} // namespace inplace

// C++17 nodiscard
template<class Container, class Predicate>
Container copy_if(const Container& container, Predicate&& predicate) {
  Container returnContainer;

  for(const auto& value : container) {
    if(invoke(predicate, value)) {
      addToContainer(returnContainer, value);
    }
  }

  return returnContainer;
}

// C++17 nodiscard
template<class Container>
Container sort(Container container) {
  inplace::sort(container);
  return container;
}

// C++17 nodiscard
template<class Container, typename Comparator>
Container sort(Container container, Comparator&& comparator) {
  inplace::sort(container, std::forward<Comparator>(comparator));
  return container;
}

// C++17 nodiscard
template<class Container, class UnaryFunction>
Container transform(Container container, UnaryFunction&& function) {
  inplace::transform(container, std::forward<UnaryFunction>(function));
  return container;
}

template<class Container, typename T>
auto find(const Container& container, const T& needle) {
  return std::find(
    std::begin(container),
    std::end(container),
    needle
  );
}

template<class Container, typename UnaryPredicate>
auto find_if(const Container& container, UnaryPredicate&& predicate) {
  return std::find_if(
    std::begin(container),
    std::end(container),
    std::forward<UnaryPredicate>(predicate)
  );
}

template<class Container>
void reverse(Container container) {
  inplace::reverse(container);
  return container;
}

template<class Container, typename T>
void remove(
  Container container,
  const T& value
) {
  inplace::remove(container, value);
  return container;
}

template<class Container, class UnaryFunction>
void remove_if(
  Container& container,
  UnaryFunction&& predicate
) {
  inplace::remove_if(container, std::forward<UnaryFunction>(predicate));
  return container;
}

template<typename UnsignedType>
std::vector<UnsignedType> iota(UnsignedType upperBound) {
  std::vector<UnsignedType> values (upperBound);

  std::iota(
    std::begin(values),
    std::end(values),
    0
  );

  return values;
}

template<class Container>
auto makeContainsPredicate(const Container& container) {
  using T = traits::getValueType<Container>;

  return [&container](const T& element) -> bool {
    return std::find(
      std::begin(container),
      std::end(container),
      element
    ) != std::end(container);
  };
}

} // namespace temple

#endif
