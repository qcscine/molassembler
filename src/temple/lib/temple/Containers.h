#ifndef INCLUDE_TEMPLATE_MAGIC_CONTAINERS_H
#define INCLUDE_TEMPLATE_MAGIC_CONTAINERS_H

#include "AddToContainer.h"
#include "constexpr/TupleType.h"

#include <algorithm>
#include <numeric>
#include <set>
#include <map>
#include <string>

/*! @file
 *
 * Provides a slew of pseudo-functional-style composable functions to ease
 * manipulation of containers with lambdas. Functions are typically geared
 * towards the minimum set of requirements expected of the containers in order
 * to function well. Most functions will work well with many STL containers,
 * and custom containers that fulfill the function's explicit requirements.
 */

namespace temple {

/* Header */

template<typename Container>
void sort(Container& container);

template<typename Container, typename Comparator>
void sort(Container& container, Comparator&& comparator);

template<typename Container, typename T>
auto find(const Container& container, const T& needle);

template<typename Container, typename UnaryPredicate>
auto find_if(const Container& container, UnaryPredicate&& predicate);

template<typename T>
std::vector<T> iota(T size);

//! Composability improvement - returns the call to the size member
template<typename Container>
auto size(const Container& container);

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
);

//! Map specialization for std::array
template<
  typename T,
  long unsigned size,
  class UnaryFunction
> auto map(
  const std::array<T, size>& container,
  UnaryFunction&& function
);

// TODO these should be replaced with a range adaptor, not be present as individual functions
// perhaps implement these as a special case of the MemberFetcher
//! Map specialization for std::map that maps the keys, returns a vector
template<typename T, typename U, class UnaryFunction>
auto mapKeys(
  const std::map<T, U>& map,
  UnaryFunction&& function
);

// TODO these should be replaced with a range adaptor, not be present as individual functions
// perhaps implement these as a special case of the MemberFetcher
//! Map specialization for std::map, which maps the values, returns a vector
template<typename T, typename U, class UnaryFunction>
auto mapValues(
  const std::map<T, U>& map,
  UnaryFunction&& function
);

/*!
* Maps the contents of a Container to a vector using a unary function.
*
* Requires:
* - Container must implement begin and end forward-iterators
* - UnaryFunction must be unary and callable with the type pointed to by
*   the iterators
*/
template<
  class Container,
  class UnaryFunction
> auto mapToVector(
  const Container& container,
  UnaryFunction&& function
);

/*!
 * Composable pairwise map function. Instead of a unary function acting on
 * every element of the container, this takes pairs of successive elements.
 * The returned type is a container of the same type as the input, containing
 * elements of the type that the binary function returns.
 */
template<
  class BinaryFunction,
  typename T,
  template<typename, typename...> class Container,
  template<typename> class ... Dependents
> auto mapSequentialPairs(
  const Container<T, Dependents<T>...>& container,
  BinaryFunction&& function
);

//! Composable pairwise map function specialization for array types
template<
  class BinaryFunction,
  template<typename, std::size_t> class ArrayType,
  typename T,
  std::size_t N
> auto mapSequentialPairs(
  const ArrayType<T, N>& array,
  BinaryFunction&& function
);


/*!
 * Takes a container and maps all possible pairs of its contents into a new
 * container of the same type.
 */
template<
  class BinaryFunction,
  typename T,
  template<typename, typename...> class Container,
  template<typename> class ... Dependents
> auto mapAllPairs(
  const Container<T, Dependents<T>...>& container,
  BinaryFunction&& function
);

/*!
 * Zip mapping. Always returns a vector containing the type returned by
 * the binary function.
 */
template<
  class ContainerT,
  class ContainerU,
  class BinaryFunction
> auto zipMap(
  const ContainerT& containerT,
  const ContainerU& containerU,
  BinaryFunction&& function
);

/*!
 * Composable reduce function. Requires that container implements begin and end
 * iterators pointing to Ts. BinaryFunction must take two Ts and return a T.
 */
template<typename Container, typename T, class BinaryFunction>
std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    T
  >::value,
  T
> reduce(
  const Container& container,
  T init,
  BinaryFunction&& binaryFunction
);

//! Takes a container and calls a function on all possible pairs of its contents
template<
  class Container,
  class BinaryFunction
> void forAllPairs(
  const Container& container,
  BinaryFunction&& function
);

/*!
 * Takes two containers and calls a function on all possible cross-container
 * pairs
 */
template<
  class ContainerT,
  class ContainerU,
  class BinaryFunction
> void forAllPairs(
  const ContainerT& containerT,
  const ContainerU& containerU,
  BinaryFunction&& function
);

//! Makes std::accumulate composable for most STL containers
template<
  class BinaryFunction,
  typename ReturnType,
  typename T,
  template<typename, typename...> class Container,
  class ... Dependents
> ReturnType accumulate(
  const Container<T, Dependents...>& container,
  ReturnType&& init,
  BinaryFunction&& function
) PURITY_WEAK;

//! Special fix for std::array accumulate composability
template<
  class BinaryFunction,
  typename ReturnType,
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> ReturnType accumulate(
  const ArrayType<T, size>& array,
  ReturnType&& init,
  BinaryFunction&& function
);

//! Returns a count of some type in a container
template<typename T, class Container> unsigned count(
  const Container& container,
  const T& toCount
);

namespace detail {

template<typename Vector> Vector concatenateHelper(Vector& vector) {
  return vector;
}

template<typename Vector, typename Container, typename ... Containers> Vector concatenateHelper(
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
template<typename... Containers> auto concatenate(
  const Containers& ... containers
) {
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

//!  Cast the entire data of a container
template<
  typename U,
  typename T,
  template<typename, typename...> class Container,
  template<typename> class ... Dependents
> Container<U, Dependents<U>...> cast(
  const Container<T, Dependents<T>...>& container
);

/*!
 * Reverses a container. Requires that the container implements begin and
 * end bidirectional iterators and a copy constructor.
 */
template<typename Container>
Container reverse(const Container& container);

/*!
 * Condenses an iterable container into a comma-separated string of string
 * representations of its contents. Requires container iterators to satisfy
 * ForwardIterators and the contained type to be a valid template
 * argument for std::to_string, which in essence means this works only for (a
 * few) STL containers and (most) built-in datatypes.
 */
template<class Container>
std::enable_if_t<
  !std::is_same<
    traits::getValueType<Container>,
    std::string
  >::value,
  std::string
> condenseIterable(
  const Container& container,
  std::string joiningChar = ", "
);

//! Version for strings
template<class Container> std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    std::string
  >::value,
  std::string
> condenseIterable(
  const Container& container,
  std::string joiningChar = ", "
);

/*!
 * Split a container's values by a mapping function whose return value elements
 * are compared with. Matching mapped values are grouped and returned in a
 * ragged 2D vector.
 */
template<class Container, class UnaryFunction>
std::vector<
  std::vector<
    traits::getValueType<Container>
  >
> groupByMapping(
  const Container& container,
  UnaryFunction&& function
);

/*!
 * Split a container's values by a binary comparison function. Returns a ragged
 * 2D vector.
 *
 * @note Requires that the equality comparison is symmetric and transitive!
 */
template<class Container, class BinaryFunction>
std::vector<
  std::vector<
    traits::getValueType<Container>
  >
> groupByEquality(
  const Container& container,
  BinaryFunction&& compareEqual
);

//! Creates a copy of the container with elements passing the predicate test
template<class Container, class UnaryFunction>
Container copyIf(
  const Container& container,
  UnaryFunction&& predicate
);

/*!
 * Moves all values passing a predicate test of a container into a new one of
 * the same type
 */
template<class Container, class UnaryFunction>
Container moveIf(
  Container&& container,
  UnaryFunction&& predicate
);

/* Reduction shorthands */
//! Tests if all elements of a container evaluate true against a predicate
template<class Container, class UnaryPredicate>
bool all_of(const Container& container, UnaryPredicate&& predicate);

//! Tests if all elements of a container are true
template<class Container>
bool all_of(const Container& container) PURITY_WEAK;

//! Tests if any elements of a container evaluate true against a predicate
template<class Container, class UnaryPredicate>
bool any_of(const Container& container, UnaryPredicate&& predicate);

//! Tests if any elements of a container are true
template<class Container>
bool any_of(const Container& container);

/* In-place algorithm shortcuts */
template<
  typename T,
  template<typename, typename> class Container
> void inplaceRemove(
  Container<T, std::allocator<T>>& container,
  const T& value
);

template<class Container, typename UnaryPredicate>
void inplaceRemoveIf(
  Container& container,
  UnaryPredicate&& predicate
);

/* Predicate HOFs */
/*!
 * HOF returning a predicate testing the container for whether it contains its
 * argument.
 */
template<class Container>
auto makeContainsPredicate(const Container& container);

/* Algorithms for special STL types */
/*!
 * Inverts the map. Returns a map that maps the opposite way. Be warned that
 * this will lead to loss of information if the original map has duplicate
 * mapped values.
 */
template<typename T, typename U>
std::map<U, T> invertMap(const std::map<T, U>& map);

//!  Composable set intersection
template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> setIntersection(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
);

//! Composable set union
template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> setUnion(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
);

//! Composable (symmetric) set difference
template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> setDifference(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
);

template<
  typename Container,
  typename ComparisonFunction,
  typename MappingFunction
> typename Container::const_iterator select(
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
  class NAryFunction,
  template<typename, std::size_t> class ArrayType,
  typename T,
  std::size_t N
> auto unpackArrayToFunction(
  const ArrayType<T, N>& array,
  NAryFunction&& function
);

/* Implementation ------------------------------------------------------------*/
template<typename Container>
void sort(Container& container) {
  std::sort(
    std::begin(container),
    std::end(container)
  );
}

template<typename Container, typename Comparator>
void sort(Container& container, Comparator&& comparator) {
  std::sort(
    std::begin(container),
    std::end(container),
    std::forward<Comparator>(comparator)
  );
}

template<typename Container, typename T>
auto find(const Container& container, const T& needle) {
  return std::find(
    std::begin(container),
    std::end(container),
    needle
  );
}

template<typename Container, typename UnaryPredicate>
auto find_if(const Container& container, UnaryPredicate&& predicate) {
  return std::find_if(
    std::begin(container),
    std::end(container),
    predicate
  );
}

template<typename T>
std::vector<T> iota(T size) {
  std::vector<T> a;
  a.reserve(size);

  for(T i = 0; i < size; ++i) {
    a.push_back(i);
  }

  return a;
}

template<typename Container>
auto size(
  const Container& container
) {
  return container.size();
}

template<
  class UnaryFunction,
  typename T,
  template<typename, typename...> class Container,
  template<typename> class ... Dependents
> auto map(
  const Container<T, Dependents<T>...>& container,
  UnaryFunction&& function
) {
  using U = traits::functionReturnType<UnaryFunction, T>;

  Container<U, Dependents<U>...> returnContainer;
  reserveIfPossible(returnContainer, container);

  for(const auto& element : container) {
    addToContainer(returnContainer, function(element));
  }

  return returnContainer;
}

template<
  typename T,
  long unsigned size,
  class UnaryFunction
> auto map(
  const std::array<T, size>& container,
  UnaryFunction&& function
) {
  using FunctionReturnType = decltype(
    function(
      std::declval<T>()
    )
  );

  std::array<FunctionReturnType, size> result;

  for(unsigned i = 0; i < size; i++) {
    result[i] = function(container[i]);
  }

  return result;
}

template<typename T, typename U, class UnaryFunction>
auto mapKeys(
  const std::map<T, U>& map,
  UnaryFunction&& function
) {
  using R = traits::functionReturnType<UnaryFunction, T>;

  std::vector<R> result;

  for(const auto& iterPair : map) {
    result.emplace_back(
      function(iterPair.first)
    );
  }

  return result;
}

template<typename T, typename U, class UnaryFunction>
auto mapValues(
  const std::map<T, U>& map,
  UnaryFunction&& function
) {
  using R = traits::functionReturnType<UnaryFunction, U>;

  std::vector<R> result;

  for(const auto& iterPair : map) {
    result.emplace_back(
      function(iterPair.second)
    );
  }

  return result;
}

template<typename Container, typename T, class BinaryFunction>
std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    T
  >::value,
  T
> reduce(
  const Container& container,
  T init,
  BinaryFunction&& binaryFunction
) {
  for(const auto& value: container) {
    init = binaryFunction(init, value);
  }

  return init;
}

template<
  class Container,
  class UnaryFunction
> auto mapToVector(
  const Container& container,
  UnaryFunction&& function
) {
  using U = traits::functionReturnType<UnaryFunction, traits::getValueType<Container>>;

  std::vector<U> returnVector;
  reserveIfPossible(returnVector, container);

  for(const auto& element : container) {
    returnVector.emplace_back(
      function(element)
    );
  }

  return returnVector;
}

template<
  class BinaryFunction,
  typename T,
  template<typename, typename...> class Container,
  template<typename> class ... Dependents
> auto mapSequentialPairs(
  const Container<T, Dependents<T>...>& container,
  BinaryFunction&& function
) {
  using U = decltype(
    function(
      std::declval<T>(),
      std::declval<T>()
    )
  );

  Container<U, Dependents<U>...> returnContainer;

  reserveIfPossible(
    returnContainer,
    container,
    [](const auto size) {
      return size - 1;
    }
  );

  auto leftIterator = std::begin(container);
  auto rightIterator = leftIterator; ++rightIterator;

  while(rightIterator != std::end(container)) {
    addToContainer(
      returnContainer,
      function(
        *leftIterator,
        *rightIterator
      )
    );

    ++leftIterator;
    ++rightIterator;
  }

  return returnContainer;
}

//! Composable pairwise map function specialization for array types
template<
  class BinaryFunction,
  template<typename, std::size_t> class ArrayType,
  typename T,
  std::size_t N
> auto mapSequentialPairs(
  const ArrayType<T, N>& array,
  BinaryFunction&& function
) {
  static_assert(N > 1, "No sequential pairs can be mapped in an array of size 1");
  using U = decltype(
    function(
      std::declval<T>(),
      std::declval<T>()
    )
  );

  ArrayType<U, N - 1> returnContainer;

  auto leftIterator = std::begin(array);
  auto rightIterator = leftIterator; ++rightIterator;
  auto insertIterator = std::begin(returnContainer);

  while(rightIterator != std::end(array)) {
    *insertIterator = function(
      *leftIterator,
      *rightIterator
    );

    ++leftIterator;
    ++rightIterator;
    ++insertIterator;
  }

  return returnContainer;
}

template<
  class BinaryFunction,
  typename T,
  template<typename, typename...> class Container,
  template<typename> class ... Dependents
> auto mapAllPairs(
  const Container<T, Dependents<T>...>& container,
  BinaryFunction&& function
) {
  using U = decltype(
    function(
      std::declval<T>(),
      std::declval<T>()
    )
  );

  Container<U, Dependents<U>...> returnContainer;
  reserveIfPossible(
    returnContainer,
    container,
    [](const auto size) {
      if(size < 2) {
        return 0;
      }

      return size * (size - 1) / 2;
    }
  );

  for(auto i = std::begin(container); i != std::end(container); ++i) {
    auto j = i; ++j;
    for(/* init before */; j != std::end(container); ++j) {
      addToContainer(
        returnContainer,
        function(
          *i,
          *j
        )
      );
    }
  }

  return returnContainer;
}

template<
  class ContainerT,
  class ContainerU,
  class BinaryFunction
> auto zipMap(
  const ContainerT& containerT,
  const ContainerU& containerU,
  BinaryFunction&& function
) {
  using T = traits::getValueType<ContainerT>;
  using U = traits::getValueType<ContainerU>;

  using FunctionReturnType = decltype(
    function(
      std::declval<T>(),
      std::declval<U>()
    )
  );

  std::vector<FunctionReturnType> data;

  const auto tEnd = std::end(containerT);
  const auto uEnd = std::end(containerU);

  auto tIter = std::begin(containerT);
  auto uIter = std::begin(containerU);

  while(tIter != tEnd && uIter != uEnd) {
    data.emplace_back(
      function(
        *tIter,
        *uIter
      )
    );

    ++tIter;
    ++uIter;
  }

  return data;
}

template<
  class BinaryFunction,
  typename ReturnType,
  typename T,
  template<typename, typename...> class Container,
  class ... Dependents
> ReturnType accumulate(
  const Container<T, Dependents...>& container,
  ReturnType&& init,
  BinaryFunction&& function
) {
  return std::accumulate(
    std::begin(container),
    std::end(container),
    init,
    std::forward<BinaryFunction>(function)
  );
}

template<
  class BinaryFunction,
  typename ReturnType,
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> ReturnType accumulate(
  const ArrayType<T, size>& array,
  ReturnType&& init,
  BinaryFunction&& function
) {
  return std::accumulate(
    std::begin(array),
    std::end(array),
    init,
    std::forward<BinaryFunction>(function)
  );
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

template<typename T, typename U>
std::map<U, T> invertMap(const std::map<T, U>& map) {
  std::map<U, T> flipped;

  for(const auto& mapPair : map) {
    flipped[mapPair.second] = mapPair.first;
  }

  return flipped;
}

template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> setIntersection(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
) {
  std::set<T, Comparator<T>, Allocator<T>> returnSet;

  std::set_intersection(
    std::begin(a),
    std::end(a),
    std::begin(b),
    std::end(b),
    std::inserter(returnSet, std::end(returnSet)),
    Comparator<T>()
  );

  return returnSet;
}

template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> setUnion(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
) {
  std::set<T, Comparator<T>, Allocator<T>> returnSet;

  std::set_union(
    std::begin(a),
    std::end(a),
    std::begin(b),
    std::end(b),
    std::inserter(returnSet, std::end(returnSet)),
    Comparator<T>()
  );

  return returnSet;
}

template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
>
std::set<T, Comparator<T>, Allocator<T>> setDifference(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
) {
  std::set<T, Comparator<T>, Allocator<T>> returnSet;

  std::set_symmetric_difference(
    std::begin(a),
    std::end(a),
    std::begin(b),
    std::end(b),
    std::inserter(returnSet, std::end(returnSet)),
    Comparator<T>()
  );

  return returnSet;
}

template<typename T, class Container> unsigned count(
  const Container& container,
  const T& toCount
) {
  unsigned count = 0;

  for(const auto& element : container) {
    if(element == toCount) {
      count += 1;
    }
  }

  return count;
}

template<
  typename U,
  typename T,
  template<typename, typename...> class Container,
  template<typename> class ... Dependents
> Container<U, Dependents<U>...> cast(
  const Container<T, Dependents<T>...>& container
) {
  Container<U, Dependents<U>...> casted;

  for(const T& value : container) {
    casted.emplace_back(
      static_cast<U>(value)
    );
  }

  return casted;
}

template<typename Container>
Container reverse(const Container& container) {
  Container copy = container;

  std::reverse(
    std::begin(copy),
    std::end(copy)
  );

  return copy;
}

template<class Container>
std::enable_if_t<
  !std::is_same<
    traits::getValueType<Container>,
    std::string
  >::value,
  std::string
> condenseIterable(
  const Container& container,
  std::string joiningChar
) {
  using namespace std::string_literals;

  std::string representation;

  for(auto it = std::begin(container); it != std::end(container); /*-*/) {
    representation += std::to_string(*it);
    if(++it != std::end(container)) {
      representation += joiningChar;
    }
  }

  return representation;
}

template<class Container> std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    std::string
  >::value,
  std::string
> condenseIterable(
  const Container& container,
  std::string joiningChar
) {
  using namespace std::string_literals;

  std::string representation;

  for(auto it = std::begin(container); it != std::end(container); /*-*/) {
    representation += *it;
    if(++it != std::end(container)) {
      representation += joiningChar;
    }
  }

  return representation;
}

template<class Container, class UnaryFunction>
std::vector<
  std::vector<
    traits::getValueType<Container>
  >
> groupByMapping(
  const Container& container,
  UnaryFunction&& function
) {
  using T = traits::getValueType<Container>;
  using R = traits::functionReturnType<UnaryFunction, T>;

  std::vector<
    std::vector<T>
  > groups;

  std::map<R, unsigned> indexMap;

  for(auto iter = std::begin(container); iter != std::end(container); ++iter) {
    auto ret = function(*iter);
    if(indexMap.count(ret) == 0) {
      indexMap[ret] = groups.size();

      groups.emplace_back(
        std::vector<T> {*iter}
      );
    } else {
      groups.at(
        indexMap.at(ret)
      ).push_back(*iter);
    }
  }

  return groups;
}

template<class Container, class BinaryFunction>
std::vector<
  std::vector<
    traits::getValueType<Container>
  >
> groupByEquality(
  const Container& container,
  BinaryFunction&& compareEqual
) {
  using T = traits::getValueType<Container>;

  std::vector<
    std::vector<T>
  > groups;

  for(auto iter = std::begin(container); iter != std::end(container); ++iter) {
    bool foundEqual = false;
    for(auto& group : groups) {

      if(compareEqual(*iter, *std::begin(group))) {
        group.push_back(*iter);
        foundEqual = true;
        break;
      }
    }

    if(!foundEqual) {
      groups.emplace_back(
        std::vector<T> {*iter}
      );
    }
  }

  return groups;
}

template<class Container, class UnaryFunction>
Container copyIf(
  const Container& container,
  UnaryFunction&& predicate
) {
  Container returnContainer;
  reserveIfPossible(
    returnContainer,
    container
  );

  for(const auto& elem : container) {
    if(predicate(elem)) {
      addToContainer(returnContainer, elem);
    }
  }

  return returnContainer;
}

template<class Container, class UnaryFunction>
Container moveIf(
  Container&& container,
  UnaryFunction&& predicate
) {
  Container returnContainer;
  reserveIfPossible(returnContainer, container);

  for(auto& elem : container) {
    if(predicate(elem)) {
      addToContainer(returnContainer, std::move(elem));
    }
  }

  return returnContainer;
}

template<class Container, class UnaryPredicate>
bool all_of(const Container& container, UnaryPredicate&& predicate) {
  for(const auto& element : container) {
    if(!predicate(element)) {
      return false;
    }
  }

  return true;
}

template<class Container>
bool all_of(const Container& container) {
  return accumulate(
    container,
    true,
    std::logical_and<bool>()
  );
}

template<class Container, class UnaryPredicate>
bool any_of(const Container& container, UnaryPredicate&& predicate) {
  for(const auto& element : container) {
    if(predicate(element)) {
      return true;
    }
  }

  return false;
}

template<class Container>
bool any_of(const Container& container) {
  return accumulate(
    container,
    false,
    std::logical_or<bool>()
  );
}

template<
  typename T,
  template<typename, typename> class Container
> void inplaceRemove(
  Container<T, std::allocator<T>>& container,
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

template<class Container, typename UnaryPredicate>
void inplaceRemoveIf(
  Container& container,
  UnaryPredicate&& predicate
) {
  container.erase(
    std::remove_if(
      std::begin(container),
      std::end(container),
      predicate
    ),
    std::end(container)
  );
}

template<
  class Container,
  class BinaryFunction
> void forAllPairs(
  const Container& container,
  BinaryFunction&& function
) {
  for(auto i = std::begin(container); i != std::end(container); ++i) {
    auto j = i; ++j;
    for(/* init before */; j != std::end(container); ++j) {
      function(
        *i,
        *j
      );
    }
  }
}

template<
  class ContainerT,
  class ContainerU,
  class BinaryFunction
> void forAllPairs(
  const ContainerT& containerT,
  const ContainerU& containerU,
  BinaryFunction&& function
) {
  for(auto i = std::begin(containerT); i != std::end(containerT); ++i) {
    for(auto j = std::begin(containerU); j != std::end(containerU); ++j) {
      function(*i, *j);
    }
  }
}

namespace detail {

template<
  class NAryFunction,
  template<typename, std::size_t> class ArrayType,
  typename T,
  std::size_t N,
  std::size_t ... Inds
> auto unpackArrayToFunctionHelper(
  const ArrayType<T, N>& array,
  NAryFunction&& function,
  std::index_sequence<Inds...>
) {
  return function(
    array.at(Inds)...
  );
}


} // namespace detail

template<
  class NAryFunction,
  template<typename, std::size_t> class ArrayType,
  typename T,
  std::size_t N
> auto unpackArrayToFunction(
  const ArrayType<T, N>& array,
  NAryFunction&& function
) {
  return detail::unpackArrayToFunctionHelper(
    array,
    std::forward<NAryFunction>(function),
    std::make_index_sequence<N>{}
  );
}

} // namespace temple

#endif
