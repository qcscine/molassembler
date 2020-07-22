/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Functional-style container-related algorithms
 *
 * Supplies a large number of functional-style algorithm definitions. Most are
 * STL algorithm shortcuts but many are composable with the range adaptors
 * provided in Adaptors/ using the function call template functions in
 * Invoke.h.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_FUNCTIONAL_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_FUNCTIONAL_H

#include "Molassembler/Temple/Invoke.h"
#include "Molassembler/Temple/AddToContainer.h"

#include <numeric>

namespace Scine {
namespace Molassembler {
namespace Temple {

template<
  class Container,
  std::enable_if_t<Traits::hasSize<Container>::value, int> = 0
> auto size(const Container& container) {
  return container.size();
}

template<
  class Container,
  std::enable_if_t<!Traits::hasSize<Container>::value, int> = 0
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
 * @complexity{@math{\Theta(N)}}
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

//! @overload
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

//! @overload
template<
  template<typename, typename> class PairType,
  typename T,
  class UnaryFunction
> auto map_stl(
  const PairType<T, T>& pair,
  UnaryFunction&& function
) {
  using U = decltype(
    invoke(function, pair.first)
  );

  return PairType<U, U> {
    function(pair.first),
    function(pair.second)
  };
}

/*! @brief Maps all values of a container using a unary function. Returns a vector
 *
 * @complexity{@math{\Theta(N)}}
 */
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
    returnContainer.push_back(
      invoke(function, value)
    );
  }

  return returnContainer;
}

template<class Pairlike, typename Unary>
auto mapHomogeneousPairlike(Pairlike p, Unary&& f) {
  return std::make_pair(f(p.first), f(p.second));
}

//! Apply a callable for all values of a container
template<class Container, class Callable>
void forEach(
  const Container& container,
  Callable&& callable
) {
  for(const auto& value : container) {
    invoke(callable, value);
  }
}

/**
 * @brief Select a value from a container
 *
 * @tparam Container Container type
 * @tparam ComparisonFunction Comparison between members of the value
 * @tparam MappingFunction Map function from ValueType to ComparableType
 * @param container container with values
 * @param comparator Comparator instance
 * @param mapFunction Map function
 *
 * @returns An iterator to the value whose mapped value is smallest value
 * according to the comparator
 */
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

/**
 * @brief Accumulate shorthand
 *
 * @complexity{@math{\Theta(N)}}
 */
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

/*! @brief all_of shorthand
 *
 * @complexity{@math{O(N)}}
 */
template<class Container, class UnaryPredicate = Functor::Identity>
bool all_of(const Container& container, UnaryPredicate&& predicate = UnaryPredicate {}) {
  for(const auto& element : container) {
    if(!invoke(predicate, element)) {
      return false;
    }
  }

  return true;
}

/*! @brief any_of shorthand
 *
 * @complexity{@math{O(N)}}
 */
template<class Container, class UnaryPredicate = Functor::Identity>
bool any_of(const Container& container, UnaryPredicate&& predicate = UnaryPredicate {}) {
  for(const auto& element : container) {
    if(invoke(predicate, element)) {
      return true;
    }
  }

  return false;
}

//! Calls std::sort on a container
template<class Container>
void sort(Container& container) {
  std::sort(
    std::begin(container),
    std::end(container)
  );
}

//! Calls std::sort with a custom comparator on a container
template<class Container, typename Comparator>
void sort(Container& container, Comparator&& comparator) {
  std::sort(
    std::begin(container),
    std::end(container),
    std::forward<Comparator>(comparator)
  );
}

//! In-place reversal
template<class Container>
void reverse(Container& container) {
  std::reverse(
    std::begin(container),
    std::end(container)
  );
}

//! Applies erase-remove idiom in-place to container
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

//! Applies erase-remove idiom in-place to container with a predicate function
template<class Container, class UnaryFunction>
void remove_if(
  Container& container,
  UnaryFunction&& predicate
) {
  container.erase(
    std::remove_if(
      std::begin(container),
      std::end(container),
      std::forward<UnaryFunction>(predicate)
    ),
    std::end(container)
  );
}

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
Container sorted(Container container) {
  sort(container);
  return container;
}

// C++17 nodiscard
template<class Container, typename Comparator>
Container sorted(Container container, Comparator&& comparator) {
  sort(container, std::forward<Comparator>(comparator));
  return container;
}

//! @brief std::find shorthand
template<class Container, typename T>
auto find(const Container& container, const T& needle) {
  return std::find(
    std::begin(container),
    std::end(container),
    needle
  );
}

//! @brief std::find_if shorthand
template<class Container, typename UnaryPredicate>
auto find_if(const Container& container, UnaryPredicate&& predicate) {
  return std::find_if(
    std::begin(container),
    std::end(container),
    std::forward<UnaryPredicate>(predicate)
  );
}

//! @brief vector iota shorthand
template<typename T>
std::vector<T> iota(unsigned upperBound) {
  std::vector<T> values(upperBound);

  std::iota(
    std::begin(values),
    std::end(values),
    T(0)
  );

  return values;
}

//! @brief Creates a predicate that uses std::find for the container's value type
template<class Container>
auto makeContainsPredicate(const Container& container) {
  using T = Traits::getValueType<Container>;

  return [&container](const T& element) -> bool {
    return std::find(
      std::begin(container),
      std::end(container),
      element
    ) != std::end(container);
  };
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
