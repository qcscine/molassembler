#ifndef TEMPLATE_MAGIC_INCLUDE_H
#define TEMPLATE_MAGIC_INCLUDE_H

#include <functional>
#include <algorithm>
#include <set>
#include "boost/optional.hpp"

namespace TemplateMagic {

/*! 
 * To help with creating consistent logical operators for multi-component
 * structs / classes. E.g.::
 *   
 *   struct Foo {
 *     unsigned a, b, c;
 *
 *     bool operator < (const Foo& other) {
 *       return componentSmaller(a, other.a).value_or( // sort by a first
 *         componentSmaller(b, other.b).value_or( // then by b
 *           componentSmaller(c, other.c).value_or( // then by c
 *             false
 *           )
 *         )
 *       );
 *     }
 *
 *     // the alternative can quickly become a lengthy chain and does not 
 *     // terminate as quickly if x > other.x
 *     bool operator < (const Foo& other) {
 *       return (
 *         a < other.a
 *         || (
 *           a == other.a
 *           && b < other.b
 *         ) || (
 *           a == other.a
 *           && b == other.b
 *           && c < other.c
 *         )
 *       );
 *     }
 *   };
 *
 */
template<typename T>
boost::optional<bool> componentSmaller(
  const T& a,
  const T& b
) {
  if(a < b) return true;
  else if(b < a) return false;
  else return boost::none;
}

/*! 
 * Composable sum function. Returns the type the container contains, assuming
 * monadic behavior on plus operator (value_type + value_type = value_type).
 */
template<class ContainerType>
typename ContainerType::value_type sum(const ContainerType& container) {
  return std::accumulate(
    container.begin(),
    container.end(),
    0,
    std::plus<
      typename ContainerType::value_type
    >()
  );
}

/*!
 * Composable map function. Returns the same container type, containing the
 * type that the unary function returns.
 */
template<
  typename T,
  template<typename, typename>
  class Container,
  class UnaryFunction
> 
auto map(
  const Container<T, std::allocator<T>>& container,
  UnaryFunction&& function
) {
  using FunctionReturnType = decltype(
    function(
      std::declval<T>() 
    )
  );

  Container<
    FunctionReturnType,
    std::allocator<FunctionReturnType>
  > returnContainer;

  std::transform(
    container.begin(),
    container.end(),
    std::back_inserter(returnContainer),
    function
  );

  return returnContainer;
}

/*!
 * Composable pairwise map function. Instead of a unary function acting on 
 * every element of the container, this takes two successive elements. 
 * The returned type is a container of the same type as the input, containing 
 * elements of the type that the binary function returns.
 */
template<
  typename T,
  template<typename, typename>
  class Container,
  class BinaryFunction
> 
auto pairwiseMap(
  const Container<T, std::allocator<T>>& container,
  BinaryFunction&& function
) {
  using FunctionReturnType = decltype(
    function(
      std::declval<T>(),
      std::declval<T>()
    )
  );

  Container<
    FunctionReturnType,
    std::allocator<FunctionReturnType>
  > returnContainer;

  auto inputIterator = container.begin();
  const auto stopPosition = container.end() - 1;
  
  auto outputIterator = std::inserter(returnContainer, returnContainer.end());

  while(inputIterator != stopPosition) {
    outputIterator = function(
      *inputIterator,
      *(inputIterator + 1)
    );

    inputIterator++;
  }

  return returnContainer;
}

/*!
 * Composable accumulate function.
 */
template<
  typename T,
  template<typename, typename = std::allocator<T>>
  class Container,
  class BinaryFunction,
  typename ReturnType
>
auto accumulate(
  const Container<T>& container,
  ReturnType&& init,
  BinaryFunction&& function
) {
  return std::accumulate(
    container.begin(),
    container.end(),
    init,
    function
  );
}

/*!
 * Composable set union without comparator
 */
template<typename T>
std::set<T> set_intersection(
  std::set<T> a,
  std::set<T> b
) {
  std::set<T> returnSet;
  std::set_intersection(
    a.begin(),
    a.end(),
    b.begin(),
    b.end(),
    std::inserter(returnSet, returnSet.end())
  );
  return returnSet;
}

/*!
 * Composable set union with comparator
 */
template<typename T, typename Comparator>
std::set<T, Comparator> set_intersection(
  std::set<T, Comparator> a,
  std::set<T, Comparator> b
) {
  std::set<T, Comparator> returnSet;
  std::set_intersection(
    a.begin(),
    a.end(),
    b.begin(),
    b.end(),
    std::inserter(returnSet, returnSet.end()),
    Comparator()
  );
  return returnSet;
}
/*!
 * Composable set union without comparator
 */
template<typename T>
std::set<T > set_union(
  std::set<T> a,
  std::set<T> b
) {
  std::set<T> returnSet;
  std::set_union(
    a.begin(),
    a.end(),
    b.begin(),
    b.end(),
    std::inserter(returnSet, returnSet.end())
  );
  return returnSet;
}

/*!
 * Composable set union with comparator
 */
template<typename T, typename Comparator>
std::set<T, Comparator> set_union(
  std::set<T, Comparator> a,
  std::set<T, Comparator> b
) {
  std::set<T, Comparator> returnSet;
  std::set_union(
    a.begin(),
    a.end(),
    b.begin(),
    b.end(),
    std::inserter(returnSet, returnSet.end()),
    Comparator()
  );
  return returnSet;
}

/*!
 * Identical pair map
 */
template<typename T, class UnaryFunction>
auto pair_map(
  const std::pair<T, T>& pair,
  UnaryFunction&& function
) {
  using FunctionReturnType = decltype(UnaryFunction(std::declval<T>()));
  return std::pair<FunctionReturnType, FunctionReturnType>({
    function(pair.first),
    function(pair.second)
  });
}

/*!
 * Identical pair which selector
 */
template<typename T, class UnaryPredicate>
T& which(
  const std::pair<T, T>& pair,
  UnaryPredicate&& function
) {
  auto boolPair = pair_map(
    pair,
    std::forward<UnaryPredicate>(function)
  );
  if(boolPair.first && boolPair.second) {
    throw(
      std::logic_error("Both elements of the pair fulfill the predicate!")
    );
  } else if(function(pair.first)) {
    return pair.first;
  } else if(function(pair.second)) {
    return pair.second;
  } else {
    throw(
      std::logic_error("Neither element of the pair fulfills the predicate!")
    );
  }
}

// n-argument boolean XOR with the definition: only one argument may be true
namespace detail {

  constexpr unsigned TMPSum() {
    return 0;
  }

  template<typename T1, typename... T>
  constexpr unsigned TMPSum(T1 a, T ... pack) {
    return a + TMPSum(pack ...);
  }

}

template<typename ... Bools>
constexpr bool XOR(Bools ... bools) {
  return detail::TMPSum(bools ...) == 1;
}
  
/*!
 * Composable size
 */
template<typename Container>
unsigned size(
  Container container
) {
  return container.size();
}

/* TODO not sure about the value in the following two, they are very 
 * specifically geared towards vector<bool> and have little other use
 */
template<class Container>
bool all_of(
  const Container& container
) {
  return accumulate(
    container,
    true,
    std::logical_and<bool>()
  );
}

template<class Container>
bool any_of(
  const Container& container
) {
  return accumulate(
    container,
    false,
    std::logical_or<bool>()
  );
}

/*!
 * HOF returning a predicate testing the container for whether it contains its
 * argument.
 */
template<
  typename T,
  template<typename, typename = std::allocator<T>>
  class Container
>
auto makeContainsPredicate(
  const Container<T>& container
) {
  return [&container](const T& element) -> bool {
    return std::find(
      container.begin(),
      container.end(),
      element
    ) != container.end();
  };
}

}

#endif
