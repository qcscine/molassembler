#ifndef INCLUDE_CONSTEXPR_MAGIC_CONTAINERS_H
#define INCLUDE_CONSTEXPR_MAGIC_CONTAINERS_H

#include <array>
#include <cassert>
#include <limits>

namespace ConstexprMagic {

namespace traits {

template<class Function, typename ...Args>
using functionReturnType = std::result_of_t<Function(Args...)>;

} // namespace traits

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

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size,
  class UnaryFunction
> constexpr auto map(
  const ArrayType<T, size>& array,
  UnaryFunction&& function
) {
  return mapImpl(
    array,
    std::forward<UnaryFunction>(function),
    std::make_index_sequence<size>{}
  );
}

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

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t N1,
  size_t... AIndices,
  size_t N2,
  size_t... BIndices
>
constexpr ArrayType<T, N1+N2> arrayContatenateImpl(
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
  return arrayConcatenateImpl(
    a,
    b,
    std::make_index_sequence<N1>{},
    std::make_index_sequence<N2>{}
  );
}

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

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr ArrayType<T, (size + 1)> arrayPush(
  const ArrayType<T, size>& array,
  const T& element
) {
  return arrayPushImpl(
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

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr std::enable_if_t<
  (size > 0),
  ArrayType<T, (size - 1)>
> arrayPop(const ArrayType<T, size>& array) {
  static_assert(size != 0, "arrayPop target array is already empty");

  return arrayPopImpl(
    array,
    std::make_index_sequence<size - 1>{}
  );
}


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

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t sizeA,
  size_t sizeB
> constexpr std::enable_if_t<
  sizeA < sizeB,
  bool
> arraysLess(
  const ArrayType<T, sizeA>& a __attribute__((unused)),
  const ArrayType<T, sizeB>& b __attribute__((unused))
) {
  return true;
}

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t sizeA,
  size_t sizeB
> constexpr std::enable_if_t<
  (sizeA > sizeB),
  bool
> arraysLess(
  const ArrayType<T, sizeA>& a __attribute__((unused)),
  const ArrayType<T, sizeB>& b __attribute__((unused))
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

template<
  typename T,
  class Comparator,
  template<typename, size_t> class ArrayType,
  size_t size
> constexpr size_t lowerBound(
  const ArrayType<T, size>& array,
  const T& item,
  Comparator comparator
) {
  size_t it = 0, count = array.size(), step = 0, bound = 0;

  while(count > 0) {
    it = bound;
    step = count / 2;
    it += step;
    
    if(comparator(array.at(it), item)) {
      it += 1;
      bound = it;
      count -= step + 1;
    } else {
      count = step;
    }
  }

  return bound;
}

template<
  typename T,
  template<typename, size_t> class ArrayType,
  size_t size
> constexpr std::enable_if_t<
  (size == 0),
  ArrayType<T, size + 1>
> insertIntoSorted(
  const ArrayType<T, size>& array __attribute__((unused)),
  const T& item
) {
  return {item};
}

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

template<
  typename T,
  class Comparator,
  template<typename, size_t> class ArrayType,
  size_t size
> constexpr std::enable_if_t<
  (size == 0),
  ArrayType<T, size + 1>
> insertIntoSorted(
  const ArrayType<T, size>& array __attribute__((unused)),
  const T& item,
  Comparator compare __attribute__((unused))
) {
  return {item};
}

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

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr void inPlaceReverse(
  ArrayType<T, size>& data,
  const unsigned& indexFrom,
  const unsigned& indexTo
) {
  size_t a = indexFrom, b = indexTo;
  while(a != b && a != --b) {
    inPlaceSwap(data, a++, b);
  }
}

// works with std::array only in C++17
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr bool inPlaceNextPermutation(ArrayType<T, size>& data) {
  size_t i = size - 1, j = 0, k = 0;

  while(true) {
    j = i;

    if(
      i != 0
      && data.at(--i) < data.at(j)
    ) {
      k = size;

      while(
        k != 0
        && !(
          data.at(i) < data.at(--k)
        )
      ) {
        continue;
      }

      inPlaceSwap(data, i, k);
      inPlaceReverse(data, j, size);
      return true;
    }

    if(i == 0) {
      inPlaceReverse(data, 0, size);
      return false;
    }
  }
}

template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t size
> constexpr bool inPlacePreviousPermutation(ArrayType<T, size>& data) {
  size_t i = size - 1, j = 0, k = 0;

  while(true) {
    j = i;

    if(
      i != 0
      && data.at(j) < data.at(--i)
    ) {
      k = size;

      while(
        k != 0
        && !(
          data.at(--k) < data.at(i)
        )
      ) {
        continue;
      }

      inPlaceSwap(data, i, k);
      inPlaceReverse(data, j, size);
      return true;
    }

    if(i == 0) {
      inPlaceReverse(data, 0, size);
      return false;
    }
  }
}

} // namespace ConstexprMagic

#endif
