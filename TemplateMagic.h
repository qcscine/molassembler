#ifndef TEMPLATE_MAGIC_INCLUDE_H
#define TEMPLATE_MAGIC_INCLUDE_H

#include <functional>
#include <algorithm>
#include <numeric>
#include <set>
#include <vector>
#include <map>
#include "boost/optional.hpp"

namespace TemplateMagic {
/* Logical operator work -----------------------------------------------------*/
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

/* Composition improvements --------------------------------------------------*/
//! Specialized map function for arrays
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

/*! Inverts the map. Returns a map that maps the opposite way. Be warned that 
 * this will lead to loss of information if the original map has duplicate
 * values.
 */
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

/*! Alternate implementation of zip mapping where the contained types are 
 * deduced instead of assuming the STL container form.
 */
template<
  class ContainerT,
  class ContainerU,
  class BinaryFunction
>
auto zipMapAlternate(
  const ContainerT& containerT,
  const ContainerU& containerU,
  BinaryFunction&& function
) {
  using T = decltype(
    *(
      containerT.begin()
    )
  );
  using U = decltype(
    *(
      containerU.begin()
    )
  );

  using FunctionReturnType = decltype(
    function(
      std::declval<T>(),
      std::declval<U>()
    )
  );

  std::vector<FunctionReturnType> data;

  const unsigned minSize = std::min(containerT.size(), containerU.size());

  auto tIter = containerT.begin();
  auto uIter = containerU.begin();

  for(unsigned read = 0; read < minSize; read++) {
    data.emplace_back(
      function(
        *tIter,
        *uIter
      )
    );

    tIter++;
    uIter++;
  }

  return data;
}

//!  Composable accumulate function.
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

//!  Composable set union without comparator
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
//!  Composable set union without comparator
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

//!  Composable set union with comparator
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

/* Pair helper functions -----------------------------------------------------*/
template<typename T, class UnaryFunction>
auto pair_map(
  const std::pair<T, T>& pair,
  UnaryFunction&& function
) __attribute__ ((deprecated));

//!  Identical pair map
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

template<typename T, class UnaryPredicate>
T& which(
  const std::pair<T, T>& pair,
  UnaryPredicate&& function
) __attribute__ ((deprecated));

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

/* numeric helper functions --------------------------------------------------*/
namespace numeric {

/*! 
 * Composable sum function. Returns the type the container contains, assuming
 * monadic behavior on operator + (value_type + value_type = value_type).
 * Container must implement begin and end members.
 */
template<class ContainerType>
typename std::remove_reference<
  decltype(
    *(
      std::declval<ContainerType>()
    ).begin()
  )
>::type
sum(const ContainerType& container) {
  using ValueType = typename std::remove_reference<decltype(*container.begin())>::type;

  return std::accumulate(
    container.begin(),
    container.end(),
    ValueType {0},
    std::plus<ValueType>()
  );
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
  assert(container.begin() != container.end());

  return static_cast<double>(
    sum(container)
  ) / container.size();
}

} // namespace numeric

} // namespace TemplateMagic

#endif
