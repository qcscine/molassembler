#ifndef INCLUDE_TEMPLATE_MAGIC_CONTAINERS_H
#define INCLUDE_TEMPLATE_MAGIC_CONTAINERS_H

#include "AddToContainer.h"
#include "constexpr_magic/TupleType.h"

#include <algorithm>
#include <numeric>
#include <set>
#include <map>

/*! @file
 *
 * Provides a slew of pseudo-functional-style composable functions to ease 
 * manipulation of containers with lambdas. Functions are typically geared 
 * towards the minimum set of requirements expected of the containers in order
 * to function well. Most functions will work well with many STL containers,
 * and custom containers that fulfill the function's explicit requirements.
 */

namespace TemplateMagic {

/* Header */
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

//! Map specialization for std::map that maps the keys, returns a vector
template<typename T, typename U, class UnaryFunction>
auto mapKeys(
  const std::map<T, U>& map,
  UnaryFunction&& function
);

//! Map specialization for std::map, which maps the values, returns a vector
template<typename T, typename U, class UnaryFunction>
auto mapValues(
  const std::map<T, U>& map,
  UnaryFunction&& function
);

/*!
* Maps the contents of a Container to a vector using a unary function. This
* is the variant for if the container does not implement a size member.
*
* Requires:
* - Container must implement begin and end forward-iterators
* - UnaryFunction must be unary and callable with the type pointed to by 
*   the iterators
*/
template<
  class Container,
  class UnaryFunction
> std::enable_if_t<
  !traits::hasSize<Container>::value,
  std::vector<
    traits::functionReturnType<
      UnaryFunction,
      traits::getValueType<Container>
    >
  >
> mapToVector(
  const Container& container,
  UnaryFunction&& function
);

/*!
* Maps the contents of a Container to a vector using a unary function. This
* is the variant for if the container implements the size member.
*
* Requires:
* - Container must implement begin and end forward-iterators
* - UnaryFunction must be unary and callable with the type pointed to by 
*   the iterators
*/
template<
  class Container,
  class UnaryFunction
> std::enable_if_t<
  traits::hasSize<Container>::value,
  std::vector<
    traits::functionReturnType<
      UnaryFunction,
      traits::getValueType<Container>
    >
  >
> mapToVector(
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
);

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
    vector.end(),
    container.begin(),
    container.end()
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
    ConstexprMagic::TupleType::countType<
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

/*! Condenses an iterable container into a comma-separated string of string 
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
  const std::string& joiningChar = ", "
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
  const std::string& joiningChar = ", "
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
std::vector<
  traits::getValueType<Container>
> copyIf(
  const Container& container,
  UnaryFunction&& predicate
);

/* Reduction shorthands */
//! Tests if all elements of a container are true
template<class Container>
bool all_of(const Container& container);

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

//! Composable set difference
template<
  typename T,
  template<typename> class Comparator,
  template<typename >class Allocator
> 
std::set<T, Comparator<T>, Allocator<T>> setDifference(
  const std::set<T, Comparator<T>, Allocator<T>>& a,
  const std::set<T, Comparator<T>, Allocator<T>>& b
);


/* Implementation ------------------------------------------------------------*/
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

  for(const auto& element : container) {
    TemplateMagic::addToContainer(returnContainer, function(element));
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
> std::enable_if_t<
  !traits::hasSize<Container>::value,
  std::vector<
    traits::functionReturnType<
      UnaryFunction,
      traits::getValueType<Container>
    >
  >
> mapToVector(
  const Container& container,
  UnaryFunction&& function
) {
  using U = traits::functionReturnType<UnaryFunction, traits::getValueType<Container>>;

  std::vector<U> returnVector;

  for(const auto& element : container) {
    returnVector.emplace_back(
      function(element)
    );
  }

  return returnVector;
}

template<
  class Container,
  class UnaryFunction
> std::enable_if_t<
  traits::hasSize<Container>::value,
  std::vector<
    traits::functionReturnType<
      UnaryFunction,
      traits::getValueType<Container>
    >
  >
> mapToVector(
  const Container& container,
  UnaryFunction&& function
) {
  using U = traits::functionReturnType<UnaryFunction, traits::getValueType<Container>>;

  std::vector<U> returnVector;
  returnVector.reserve(container.size());

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

  auto leftIterator = container.begin();
  auto rightIterator = leftIterator; ++rightIterator;

  while(rightIterator != container.end()) {
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

  for(auto i = container.begin(); i != container.end(); ++i) {
    auto j = i; ++j;
    for(/* init before */; j != container.end(); ++j) {
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

  const auto tEnd = containerT.end();
  const auto uEnd = containerU.end();

  auto tIter = containerT.begin();
  auto uIter = containerU.begin();

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
    container.begin(),
    container.end(),
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
    array.begin(),
    array.end(),
    init,
    std::forward<BinaryFunction>(function)
  );
}

template<class Container>
auto makeContainsPredicate(const Container& container) {
  using T = traits::getValueType<Container>;

  return [&container](const T& element) -> bool {
    return std::find(
      container.begin(),
      container.end(),
      element
    ) != container.end();
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
    a.begin(),
    a.end(),
    b.begin(),
    b.end(),
    std::inserter(returnSet, returnSet.end()),
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

  std::set_intersection(
    a.begin(),
    a.end(),
    b.begin(),
    b.end(),
    std::inserter(returnSet, returnSet.end()),
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
    a.begin(),
    a.end(),
    b.begin(),
    b.end(),
    std::inserter(returnSet, returnSet.end()),
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
    copy.begin(),
    copy.end()
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
  const std::string& joiningChar
) {
  using namespace std::string_literals;

  std::string representation;

  for(auto it = container.begin(); it != container.end(); /*-*/) {
    representation += std::to_string(*it);
    if(++it != container.end()) {
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
  const std::string& joiningChar
) {
  using namespace std::string_literals;

  std::string representation;

  for(auto it = container.begin(); it != container.end(); /*-*/) {
    representation += *it;
    if(++it != container.end()) {
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

  for(auto iter = container.begin(); iter != container.end(); ++iter) {
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

  for(auto iter = container.begin(); iter != container.end(); ++iter) {
    bool foundEqual = false;
    for(auto& group : groups) {

      if(compareEqual(*iter, *group.begin())) {
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
std::vector<
  traits::getValueType<Container>
> copyIf(
  const Container& container,
  UnaryFunction&& predicate
) {
  std::vector<
    traits::getValueType<Container>
  > ret;

  for(const auto& elem : container) {
    if(predicate(elem)) {
      addToContainer(ret, elem);
    }
  }

  return ret;
}

template<class Container>
bool all_of(const Container& container) {
  return accumulate(
    container,
    true,
    std::logical_and<bool>()
  );
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
      container.begin(),
      container.end(),
      value
    ),
    container.end()
  );
}

template<class Container, typename UnaryPredicate>
void inplaceRemoveIf(
  Container& container,
  UnaryPredicate&& predicate
) {
  container.erase(
    std::remove_if(
      container.begin(),
      container.end(),
      predicate
    ),
    container.end()
  );
}

template<
  class Container,
  class BinaryFunction
> void forAllPairs(
  const Container& container,
  BinaryFunction&& function
) {
  for(auto i = container.begin(); i != container.end(); ++i) {
    auto j = i; ++j;
    for(/* init before */; j != container.end(); ++j) {
      function(
        *i,
        *j
      );
    }
  }
}

} // namespace TemplateMagic

#endif
