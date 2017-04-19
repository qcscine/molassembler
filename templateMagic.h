#ifndef TEMPLATE_MAGIC_INCLUDE_H
#define TEMPLATE_MAGIC_INCLUDE_H

#include <functional>
#include <algorithm>
#include <numeric>
#include <set>
#include <map>
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
    typename ContainerType::value_type(0),
    std::plus<
      typename ContainerType::value_type
    >()
  );
}

template<
  typename T,
  long unsigned size,
  class UnaryFunction
>
auto map(
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

/*!
 * Composable map function. Returns the same container type, containing the
 * type that the unary function returns.
 */
template<
  typename T,
  template<typename, typename> class Container,
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

template<typename T, typename U>
std::map<U, T> invertMap(const std::map<T, U>& map) {
  std::map<U, T> flipped;
  
  for(const auto& mapPair : map) {
    flipped[mapPair.second] = mapPair.first;
  }

  return flipped;
}

/*!
 * Composable pairwise map function. Instead of a unary function acting on 
 * every element of the container, this takes pairs of successive elements. 
 * The returned type is a container of the same type as the input, containing 
 * elements of the type that the binary function returns.
 */
template<
  typename T,
  template<typename, typename> class Container,
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
 * Zipping mapper for same container primitives
 */
template<
  typename T,
  typename U,
  template<typename, typename> class Container,
  class BinaryFunction
>
auto zipMap(
  const Container<T, std::allocator<T>>& a,
  const Container<U, std::allocator<T>>& b,
  BinaryFunction&& function
) {
  using FunctionReturnType = decltype(
    function(
      std::declval<T>(),
      std::declval<T>()
    )
  );

  Container<FunctionReturnType, std::allocator<FunctionReturnType>> data;

  const unsigned minSize = std::min(a.size(), b.size());

  for(unsigned i = 0; i < minSize; i++) {
    data.push_back(
      function(
        a[i],
        b[i]
      )
    );
  }

  return data;
}

/*!
 * Composable accumulate function.
 */
template<
  typename T,
  template<typename, typename = std::allocator<T>> class Container,
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

//! Identical pair which selector
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

//! Template parameter-pack exclusive or of booleans
template<typename ... Bools>
constexpr bool XOR(Bools ... bools) {
  return detail::TMPSum(bools ...) == 1;
}
  
//!  Composable size function
template<typename Container>
unsigned size(
  Container container
) {
  return container.size();
}

//! Tests if all elements of a container are true
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

//! Tests if any elements of a container are true
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
template<class Container>
auto makeContainsPredicate(
  const Container& container
) {
  using T = decltype(*container.begin());

  return [&container](const T& element) -> bool {
    return std::find(
      container.begin(),
      container.end(),
      element
    ) != container.end();
  };
}

/*! Takes a container and maps all possible pairs of its contents into a new 
 * container of the same type.
 */
template<
  typename T,
  template<typename, typename> class Container,
  class BinaryFunction
> 
auto allPairsMap(
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

  auto outputIterator = std::inserter(returnContainer, returnContainer.end());

  for(auto i = container.begin(); i != container.end(); i++) {
    for(auto j = i + 1; j != container.end(); j++) {
      outputIterator = function(
        *i,
        *j
      );
    }
  }

  return returnContainer;
}

//! Returns a count of some type in a container
template<
  typename T,
  template<typename, typename> class Container
> unsigned count(
  const Container<T, std::allocator<T>>& container,
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

/*! Condenses an iterable container into a comma-separated string of string 
 * representations of its contents. Requires container iterators to satisfy
 * BidirectionalIterators and the contained type to be a valid template
 * argument for std::to_string, which in essence means this works only for (a
 * few) STL containers and (most) built-in datatypes.
 */
template<
  typename T,
  template<typename, typename> class Container
> std::string condenseIterable(
  const Container<T, std::allocator<T>>& container
) {
  using namespace std::string_literals;

  std::string representation;

  for(auto it = container.begin(); it != container.end(); it++) {
    representation += std::to_string(*it);
    if(it != container.end() - 1) {
      representation += ", "s;
    }
  }

  return representation;
}

}

#endif
