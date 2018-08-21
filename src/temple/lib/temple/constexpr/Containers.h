#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_CONTAINERS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_CONTAINERS_H

#include "temple/Preprocessor.h"

#include <array>
#include <limits>

/*! @file
 *
 * Provides constexpr functional-style modification of container elements
 */

namespace temple {

namespace traits {

template<class Function, typename ...Args>
using functionReturnType = std::result_of_t<Function(Args...)>;

} // namespace traits

namespace detail {

//!  Implementation of mapping with a unary function for any array type.
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size,
  class UnaryFunction,
  size_t ... Inds
> constexpr auto mapImpl(
  const ArrayType<T, size>& array,
  UnaryFunction&& function,
  std::index_sequence<Inds...>
) {
  return ArrayType<
    traits::functionReturnType<UnaryFunction, T>,
    size
  > {
    function(array.at(Inds))...
  };
}

} // namespace detail

//! Maps all elements of any array-like container with a unary function
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size,
  class UnaryFunction
> constexpr auto map(
  const ArrayType<T, size>& array,
  UnaryFunction&& function
) {
  return detail::mapImpl(
    array,
    std::forward<UnaryFunction>(function),
    std::make_index_sequence<size>{}
  );
}

/*!
 * Reduction of any array-like container with a binary function. This function
 * is somewhat restricted, in that the type of the reduction must be identical
 * to the value type of the array-like container. The binary function must
 * have the signature T(const T& reduction, const T& element).
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size,
  class BinaryFunction
> constexpr T reduce(
  const ArrayType<T, size>& array,
  T init,
  BinaryFunction&& reduction
) {
  for(const auto& element : array) {
    init = reduction(init, element);
  }

  return init;
}

/*!
 * Summation of all elements of an array-like class. Requires operator + and
 * zero-initialization of the contained type
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr T sum(const ArrayType<T, size>& array) {
  T sum {0};

  for(unsigned i = 0; i < size; ++i) {
    sum += array.at(i);
  }

  return sum;
}

namespace detail {


template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size,
  size_t ... Inds
> constexpr ArrayType<T, size> iotaHelper(
  std::index_sequence<Inds...>
) {
  return ArrayType<T, size> {
    Inds...
  };
}

} // namespace detail

//! Iota for any array type
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr std::enable_if_t<
  std::is_arithmetic<T>::value,
  ArrayType<T, size>
> iota() PURITY_STRONG;

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr std::enable_if_t<
  std::is_arithmetic<T>::value,
  ArrayType<T, size>
> iota() {
  return detail::iotaHelper<ArrayType, T, size>(
    std::make_index_sequence<size>()
  );
}

namespace detail {

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t begin,
  size_t end,
  size_t ... Inds
> constexpr ArrayType<T, (end - begin)> rangeHelper(
  std::index_sequence<Inds...>
) {
  return ArrayType<T, (end - begin)> {
    (begin + Inds)...
  };
}

} // namespace detail

//! Range for any array type
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t begin, // inclusive
  size_t end // exclusive
> constexpr ArrayType<T, (end - begin)> range() {
  return detail::rangeHelper<ArrayType, T, begin, end>(
    std::make_index_sequence<(end - begin)>()
  );
}

namespace detail {

//! Implementation helper of array-like type concatenation.
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t N1,
  size_t... AIndices,
  size_t N2,
  size_t... BIndices
>
constexpr ArrayType<T, N1+N2> arrayConcatenateImpl(
  const ArrayType<T, N1>& a,
  const ArrayType<T, N2>& b,
  std::index_sequence<AIndices...>,
  std::index_sequence<BIndices...>
) {
  return {
    a.at(AIndices)...,
    b.at(BIndices)...
  };
}

} // namespace detail

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t N1,
  size_t N2
>
constexpr ArrayType<T, N1+N2> arrayConcatenate(
  const ArrayType<T, N1>& a,
  const ArrayType<T, N2>& b
) {
  return detail::arrayConcatenateImpl(
    a,
    b,
    std::make_index_sequence<N1>{},
    std::make_index_sequence<N2>{}
  );
}

namespace detail {

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t concatenatedSize
> constexpr ArrayType<T, concatenatedSize> concatenateHelper(
  const ArrayType<T, concatenatedSize>& concatenated
) {
  return concatenated;
}

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t concatenatedSize,
  size_t curSize,
  size_t ... Ns
> constexpr auto concatenateHelper(
  const ArrayType<T, concatenatedSize>& concatenated,
  const ArrayType<T, curSize>& array,
  const ArrayType<T, Ns>& ... arrays
) {
  return concatenateHelper(
    arrayConcatenate(
      concatenated,
      array
    ),
    arrays...
  );
}

} // namespace detail

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t N,
  size_t ... Ns
> constexpr auto arrayConcatenate(
  const ArrayType<T, N>& startingArray,
  const ArrayType<T, Ns>& ... remainingArrays
) {
  return detail::concatenateHelper(startingArray, remainingArrays...);
}

namespace detail {

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size,
  size_t... Indices
> constexpr ArrayType<T, (size + 1)> arrayPushImpl(
  const ArrayType<T, size>& array,
  const T& element,
  std::index_sequence<Indices...>
) {
  return {
    array.at(Indices)...,
    element
  };
}

} // namespace detail

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr ArrayType<T, (size + 1)> arrayPush(
  const ArrayType<T, size>& array,
  const T& element
) {
  return detail::arrayPushImpl(
    array,
    element,
    std::make_index_sequence<size>{}
  );
}

// C++17
/*template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr ArrayType<T, (size + 1)> arrayPush(
  const ArrayType<T, size>& array,
  const T& element
) {
  ArrayType<T, (size + 1)> newArray {};

  for(size_t i = 0; i < size; ++i) {
    newArray.at(i) = array.at(i);
  }

  newArray.at(size) = element;
  return newArray;
}*/

namespace detail {

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size,
  size_t... ArrayIndices
> constexpr ArrayType<T, (size - 1)> arrayPopImpl(
  const ArrayType<T, size>& array,
  std::index_sequence<ArrayIndices...>
) {
  return {
    array.at(ArrayIndices)...
  };
}

} // namespace detail

/*!
 * Removes the last element from an array. Only compiles if size of array-like
 * container is greater than zero
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr std::enable_if_t<
  (size > 0),
  ArrayType<T, (size - 1)>
> arrayPop(const ArrayType<T, size>& array) {
  static_assert(size != 0, "arrayPop target array is already empty");

  return detail::arrayPopImpl(
    array,
    std::make_index_sequence<size - 1>{}
  );
}

//! Array-like container equality comparaotr
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr bool arraysEqual(
  const ArrayType<T, size>& a,
  const ArrayType<T, size>& b
) {
  for(unsigned i = 0; i < size; ++i) {
    if(a.at(i) != b.at(i)) {
      return false;
    }
  }

  return true;
}

/*!
 * Array-like container ordering comparator specialization for containers of
 * equal size
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t sizeA,
  size_t sizeB
> constexpr std::enable_if_t<
  sizeA == sizeB,
  bool
> arraysLess(
  const ArrayType<T, sizeA>& a,
  const ArrayType<T, sizeB>& b
) {
  for(unsigned i = 0; i < sizeA; ++i) {
    if(!(a.at(i) < b.at(i))) {
      return false;
    }
  }

  return true;
}

/*!
 * Array-like container ordering comparator specialization for containers of
 * unequal size, case a < b
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t sizeA,
  size_t sizeB
> constexpr std::enable_if_t<
  sizeA < sizeB,
  bool
> arraysLess(
  const ArrayType<T, sizeA>& /* a */,
  const ArrayType<T, sizeB>& /* b */
) {
  return true;
}

/*!
 * Array-like container ordering comparator specialization for containers of
 * unequal size, case a > b
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t sizeA,
  size_t sizeB
> constexpr std::enable_if_t<
  (sizeA > sizeB),
  bool
> arraysLess(
  const ArrayType<T, sizeA>& /* a */,
  const ArrayType<T, sizeB>& /* b */
) {
  return false;
}

// C++17: Replace all arraysLess functions with this single function below:
/*template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t sizeA,
  size_t sizeB
> constexpr std::enable_if_t<
  (sizeA > sizeB),
  bool
> arraysLess(
  const ArrayType<T, sizeA>& a,
  const ArrayType<T, sizeB>& b
) {
  if constexpr(sizeA < sizeB) {
    return true;
  } else if constexpr(sizeA > sizeB) {
    return false;
  } else { // sizes equal
    for(unsigned i = 0; i < sizeA; ++i) {
      if(!(a.at(i) < b.at(i))) {
        return false;
      }
    }

    return true;
  }
}*/

/*!
 * Returns an index within the passed array-like container whose element at that
 * index is not less than the passed vaue. Proceeds via binary search. 1:1
 * constexpr variant of std::lower_bound.
 */
template<
  typename T,
  class LessThanPredicate,
  class Iter
> constexpr Iter lowerBound(
  Iter bound,
  Iter last,
  const T& item,
  LessThanPredicate predicate
) {
  using DiffType = typename std::iterator_traits<Iter>::difference_type;

  Iter it = bound;
  DiffType count = last - bound, step = 0;

  while(count > 0) {
    it = bound;
    step = count / 2;
    it += step; // C++17 std::advance (differentiates between RandomAccess / linear)

    if(predicate(*it, item)) {
      ++it;
      bound = it;
      count -= step + 1;
    } else {
      count = step;
    }
  }

  return bound;
}

/*!
 * Binary searches within an ordered container. Returns an iterator to the
 * sought element if found, otherwise returns the end iterator.
 */
template<
  typename Container,
  typename T,
  class LessThanPredicate
> constexpr typename Container::const_iterator binarySearch(
  const Container& container,
  const T& item,
  LessThanPredicate predicate = std::less<T>()
) {
  auto bound = lowerBound<T, LessThanPredicate>(
    container.begin(),
    container.end(),
    item,
    predicate
  );

  /* Lower bound merely finds the first item not smaller than the sought one,
   * this does NOT mean that they're equal.
   *
   * a == b  <-->  !(a < b) && !(b < a)
   *
   * Lower bound guarantees the first condition of the right side, but we need
   * to check the second one to ensure that we have found what we seek, not just
   * merely something bigger than a.
   */
  if(predicate(item, *bound)) {
    return container.end();
  }

  return bound;
}

//! Sorted array-like container insertion specialization for an empty array
template<
  typename T,
  template<typename, size_t> class ArrayType,
  size_t size
> constexpr std::enable_if_t<
  (size == 0),
  ArrayType<T, size + 1>
> insertIntoSorted(
  const ArrayType<T, size>& /* array */,
  const T& item
) {
  return {item};
}

/*!
 * Sorted array-like container insertion specialization for an array containing
 * elements. Requires move compatibility of the contained type.
 */
template<
  typename T,
  template<typename, size_t> class ArrayType,
  size_t size
> constexpr std::enable_if_t<
  (size > 0),
  ArrayType<T, size + 1>
> insertIntoSorted(
  const ArrayType<T, size>& array,
  const T& item
) {
  auto newArray = arrayPush(array, item);
  auto newItemIt = newArray.end();
  --newItemIt;

  auto prevIter = newItemIt;

  while(newItemIt != newArray.begin()) {
    prevIter = newItemIt;
    --prevIter;

    if(item < *prevIter) {
      // Perform a swap in-place
      T intermediate = std::move(*newItemIt);
      *newItemIt = std::move(*prevIter);
      *prevIter = std::move(intermediate);
      --newItemIt;
    } else {
      break;
    }
  }

  return newArray;
}

/*!
 * Sorted array-like container insertion specialization for an empty array
 * including a custom element comparator
 */
template<
  typename T,
  class Comparator,
  template<typename, size_t> class ArrayType,
  size_t size
> constexpr std::enable_if_t<
  (size == 0),
  ArrayType<T, size + 1>
> insertIntoSorted(
  const ArrayType<T, size>& /* array */,
  const T& item,
  Comparator /* compare */
) {
  return {item};
}

/*!
 * Sorted array-like container insertion specialization for an array containing
 * elements with a custom comparator.
 */
template<
  typename T,
  class Comparator,
  template<typename, size_t> class ArrayType,
  size_t size
> constexpr std::enable_if_t<
  (size > 0),
  ArrayType<T, size + 1>
> insertIntoSorted(
  const ArrayType<T, size>& array,
  const T& item,
  Comparator compare
) {
  auto newArray = arrayPush(array, item);
  auto newItemIt = newArray.end();
  --newItemIt;

  auto prevIter = newItemIt;

  while(newItemIt != newArray.begin()) {
    prevIter = newItemIt;
    --prevIter;

    if(compare(item, *prevIter)) {
      // Perform a swap in-place
      T intermediate = std::move(*newItemIt);
      *newItemIt = std::move(*prevIter);
      *prevIter = std::move(intermediate);
      --newItemIt;
    } else {
      break;
    }
  }

  return newArray;
}

//! Index-based in-place swapping of elements in an array-like container
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr void inPlaceSwap(
  ArrayType<T, size>& data,
  const size_t& a,
  const size_t& b
) {
  T intermediate = std::move(data.at(b));
  data.at(b) = std::move(data.at(a));
  data.at(a) = std::move(intermediate);
}

//! Index-based in-place reversal of elements in an array-like container
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr void inPlaceReverse(
  ArrayType<T, size>& data,
  const unsigned indexFrom,
  const unsigned indexTo
) {
  size_t a = indexFrom, b = indexTo;
  while(a != b && a != --b) {
    inPlaceSwap(data, a++, b);
  }
}

/*!
 * In-place next permutation of elements in an array-like type. 1:1 index-based
 * variant of std::next_permutation
 * NOTE: works with std::array only in C++17
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr std::enable_if_t<
  (size > 1),
  bool
> inPlaceNextPermutation(
  ArrayType<T, size>& data,
  const size_t& first,
  const size_t& last
) {
  if(!(
    first < last
    && first < size
    && last <= size
  )) {
    throw "Call parameters to inPlaceNextPermutation make no sense!";
  }

  size_t i = last - 1, j = 0, k = 0;

  while(true) {
    j = i;

    if(
      i != 0
      && data.at(--i) < data.at(j)
    ) {
      k = last;

      while(
        k != 0
        && !(
          data.at(i) < data.at(--k)
        )
      ) {
        continue;
      }

      inPlaceSwap(data, i, k);
      inPlaceReverse(data, j, last);
      return true;
    }

    if(i == first) {
      inPlaceReverse(data, first, last);
      return false;
    }
  }
}

//! Shorthand for next permutation over entire array
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr std::enable_if_t<
  (size > 1),
  bool
> inPlaceNextPermutation(ArrayType<T, size>& data) {
  return inPlaceNextPermutation(data, 0, size);
}

/*!
 * In-place previous permutation of elements in an array-like type. 1:1
 * index-based variant of std::prev_permutation
 * NOTE: works with std::array only in C++17
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr std::enable_if_t<
  (size > 1),
  bool
> inPlacePreviousPermutation(
  ArrayType<T, size>& data,
  const size_t& first,
  const size_t& last
) {
  if(!(
    first < last
    && first < size
    && last <= size
  )) {
    throw "Call parameters to inPlaceNextPermutation make no sense!";
  }

  size_t i = last - 1, j = 0, k = 0;

  while(true) {
    j = i;

    if(
      i != 0
      && data.at(j) < data.at(--i)
    ) {
      k = last;

      while(
        k != 0
        && !(
          data.at(--k) < data.at(i)
        )
      ) {
        continue;
      }

      inPlaceSwap(data, i, k);
      inPlaceReverse(data, j, last);
      return true;
    }

    if(i == first) {
      inPlaceReverse(data, first, last);
      return false;
    }
  }
}

//! Shorthand for previous permutation over entire array
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr std::enable_if_t<
  (size > 1),
  bool
> inPlacePreviousPermutation(ArrayType<T, size>& data) {
  return inPlacePreviousPermutation(data, 0, size);
}

//! Calculates the index of permutation of a container
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr size_t permutationIndex(const ArrayType<T, size>& container) {
  size_t index = 0;
  size_t position = 2;// position 1 is paired with factor 0 and so is skipped
  size_t factor = 1;

  for(size_t p = size - 2; p != std::numeric_limits<size_t>::max(); --p) {
    size_t largerSuccessors = 0;

    for(size_t q = p + 1; q < size; ++q) {
      if(container.at(p) > container.at(q)) {
        ++largerSuccessors;
      }
    }

    index += (largerSuccessors * factor);
    factor *= position;
    ++position;
  }

  return index;
}

namespace detail {
  template<class ContainerType>
  struct getValueTypeImpl {
    using type = typename std::remove_const<
      typename std::remove_reference<
        decltype(
          *(
            std::declval<ContainerType>()
          ).begin()
        )
      >::type
    >::type;
  };
} // namespace detail

template<class ContainerType>
using getValueType = typename detail::getValueTypeImpl<ContainerType>::type;

/*!
 * Checks if a container holds strictly non-decreasing sequential values, i.e.
 * its values are partially ordered.
 *
 * Partially ordered: 1, 1, 2, 3
 * Totally ordered: 1, 2, 3 (no duplicates)
 * Neither: 2, 1, 3
 */
template<
  class ContainerType,
  class LessThanComparator = std::less<
    getValueType<ContainerType>
  >
> constexpr bool isPartiallyOrdered(
  const ContainerType& container,
  LessThanComparator comparator = LessThanComparator {}
) {
  auto leftIter = container.begin();
  if(leftIter == container.end()) {
    return true;
  }

  auto rightIter = leftIter + 1;

  while(rightIter != container.end()) {
    // equivalent to *rightIter < *leftIter
    if(comparator(*rightIter, *leftIter)) {
      return false;
    }

    ++leftIter;
    ++rightIter;
  }

  return true;
}

/*!
 * Checks if a container holds strictly increasing values, i.e. its values are
 * totally ordered.
 *
 * Partially ordered: 1, 1, 2, 3
 * Totally ordered: 1, 2, 3 (no duplicates)
 * Neither: 2, 1, 3
 */
template<
  class ContainerType,
  class LessThanComparator = std::less<
    getValueType<ContainerType>
  >
> constexpr bool isTotallyOrdered(
  const ContainerType& container,
  LessThanComparator comparator = LessThanComparator {}
) {
  auto leftIter = container.begin();
  if(leftIter == container.end()) {
    return true;
  }

  auto rightIter = leftIter + 1;


  while(rightIter != container.end()) {
    // Equivalent to !(*leftIter < *rightIter
    if(!comparator(*leftIter, *rightIter)) {
      return false;
    }

    ++leftIter;
    ++rightIter;
  }

  return true;
}

} // namespace temple

#endif
