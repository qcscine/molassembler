/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Functional-style constexpr container algorithms
 *
 * Provides constexpr functional-style modification of container elements
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_CONTAINERS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_CONTAINERS_H

#include "molassembler/Temple/Preprocessor.h"

#include <array>
#include <limits>

namespace Scine {
namespace Temple {
namespace Traits {

//! Figure out the return type of calling a function
template<class Function, typename ...Args>
using functionReturnType = std::result_of_t<Function(Args...)>;

} // namespace Traits

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
  std::index_sequence<Inds...> /* inds */
) {
  return ArrayType<
    Traits::functionReturnType<UnaryFunction, T>,
    size
  > {
    function(array.at(Inds))...
  };
}

} // namespace detail

/*! @brief Maps all elements of any array-like container with a unary function
 *
 * @complexity{@math{\Theta(N)}}
 */
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

/*! @brief Reduce an array-like container with a binary function
 *
 * Reduction of any array-like container with a binary function. This function
 * is somewhat restricted, in that the type of the reduction must be identical
 * to the value type of the array-like container. The binary function must
 * have the signature T(const T& reduction, const T& element).
 *
 * @complexity{@math{\Theta(N)}}
 * @todo rename to accumulate
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

/*! @brief Sum up all elements of an array-like class
 *
 * Summation of all elements of an array-like class. Requires operator + and
 * zero-initialization of the contained type
 *
 * @complexity{@math{\Theta(N)}}
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
  std::index_sequence<Inds...> /* inds */
) {
  return ArrayType<T, size> { static_cast<T>(Inds)...  };
}

} // namespace detail

//! Iota for any array type
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> PURITY_STRONG constexpr std::enable_if_t<
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
  std::index_sequence<Inds...> /* inds */
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
  std::index_sequence<AIndices...> /* aInds */,
  std::index_sequence<BIndices...> /* bInds */
) {
  return {
    a.at(AIndices)...,
    b.at(BIndices)...
  };
}

} // namespace detail

//! Concatenation of two instances of an array-like class
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

//! Variadic concatenation of multiple array-like class instances
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

/*! @brief Array-like container lexicographic equality comparaotr
 *
 * @complexity{@math{O(N)}}
 */
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

/*! @brief Lexicographical comparison for two instances of an array-like class
 *
 * Array-like container ordering comparator specialization for containers of
 * equal size.
 *
 * @complexity{@math{O(N)}}
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

/*! @brief Less than comparison of array-like classes for mismatched sizes
 *
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

/*! @brief Less than comparison of array-like classes for mismatched sizes
 *
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

/*! @brief Constexpr lower bound algorithm from STL
 *
 * Returns an index within the passed array-like container whose element at that
 * index is not less than the passed vaue. Proceeds via binary search. 1:1
 * constexpr variant of std::lower_bound.
 *
 * @complexity{@math{\Theta(N log N)}}
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
  using DiffType = decltype(last - bound);
  static_assert(
    std::is_signed<DiffType>::value,
    "Difference between iterators must be a signed type!"
  );

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

/*! @brief Binary search an order container
 *
 * Binary searches within an ordered container. Returns an iterator to the
 * sought element if found, otherwise returns the end iterator.
 *
 * @complexity{@math{\Theta(N log N)}}
 */
template<
  class Container,
  typename T,
  class LessThanPredicate = std::less<>
> constexpr typename Container::const_iterator binarySearch(
  const Container& container,
  const T& item,
  LessThanPredicate predicate = LessThanPredicate()
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

/*! @brief Index-based in-place swapping of elements in an array-like container
 *
 * @complexity{@math{\Theta(1)}}
 */
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

/*! @brief Index-based in-place reversal of elements in an array-like container
 *
 * @complexity{@math{\Theta(N)}}
 */
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

/*! @brief In-place next permutation
 *
 * In-place next permutation of elements in an array-like type. 1:1 index-based
 * variant of std::next_permutation
 * NOTE: works with std::array only in C++17 (missing constexpr markers)
 *
 * @complexity{@math{O(N/2)}}
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
        // continue;
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

//! @overload
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

/*! @brief In-place previous permutation
 *
 * In-place previous permutation of elements in an array-like type. 1:1
 * index-based variant of std::prev_permutation
 * NOTE: works with std::array only in C++17
 *
 * @complexity{@math{O(N/2)}}
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
        // continue;
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

//! @overload
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

/*! @brief Calculates the index of permutation of a container
 *
 * @complexity{@math{\Theta(N^2)}}
 */
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

//! Figures out the value type of a container via its iterators
template<class ContainerType>
using getValueType = typename detail::getValueTypeImpl<ContainerType>::type;

/*! @brief Checks if a container is partially ordered
 *
 * Checks if a container holds strictly non-decreasing sequential values, i.e.
 * its values are partially ordered.
 *
 * Partially ordered: 1, 1, 2, 3
 * Totally ordered: 1, 2, 3 (no duplicates)
 * Neither: 2, 1, 3
 *
 * @complexity{@math{O(N)}}
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

/*! @brief Checks if the container hold strictly increasing values
 *
 * Checks if a container holds strictly increasing values, i.e. its values are
 * totally ordered.
 *
 * Partially ordered: 1, 1, 2, 3
 * Totally ordered: 1, 2, 3 (no duplicates)
 * Neither: 2, 1, 3
 *
 * @complexity{@math{O(N)}}
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

} // namespace Temple
} // namespace Scine

#endif
