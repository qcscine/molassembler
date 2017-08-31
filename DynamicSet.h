#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_SET_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_SET_H

#include <set>

#include "Containers.h"
#include "DynamicArray.h"

/*! @file
 *
 * A constexpr fixed-maximum-size managed set so that the type signature does
 * not change upon element insertion and deletion. STL parallel is std::set, 
 * with the difference that here the maximum number of elements must be known at
 * compile time.
 */

namespace ConstexprMagic {

template<
  typename T,
  size_t nItems,
  class Comparator = std::less<T>
> class DynamicSet {
private:
  /* C++17 -> switch to std::array, reference operator[] and reference at()
   * become constexpr, remember to change constIterator type name
   */
  using ArrayType = ConstexprMagic::DynamicArray<T, nItems>;

  ArrayType _items;
  Comparator _comparator;

public:
  constexpr DynamicSet() {}

  //! Warning: These constructors expect ordered arrays!
  constexpr DynamicSet(const ArrayType& items) : _items(items) {}
  constexpr DynamicSet(ArrayType&& items) : _items(items) {}

  constexpr bool contains(const T& item) const {
    // C++17 turn into if constexpr-else
    if(nItems == 0) {
      return false;
    }

    auto bound = lowerBound<T, Comparator>(_items, item, _comparator);
    return (
      bound != _items.size()
      && !(_comparator(item, _items.at(bound)))
    );
  }

  //! Insertion can involve up to 3N moves of type T.
  /* TODO could be better, use lowerBound to find insert position, then move 
   * linearly towards end instead of triangularly for each swap
   */
  constexpr void insert(const T& item) {
    // add onto end
    _items.push_back(item);

    // re-sort _items
    auto newItemIt = _items.end();
    --newItemIt;

    auto prevIter = newItemIt;

    while(newItemIt != _items.begin()) {
      prevIter = newItemIt;
      --prevIter;

      if(_comparator(item, *prevIter)) {
        // Perform a swap in-place
        T intermediate = std::move(*newItemIt);
        *newItemIt = std::move(*prevIter);
        *prevIter = std::move(intermediate);

        --newItemIt;
      } else {
        break;
      }
    }
  }

  constexpr void clear() {
    _items.clear();
  }


  constexpr typename ArrayType::constIterator begin() const {
    return _items.begin();
  }

  constexpr typename ArrayType::constIterator end() const {
    return _items.end();
  }

  constexpr size_t size() const {
    return _items.size();
  }

  constexpr bool operator == (const DynamicSet& other) const {
    return _items == other._items;
  }

  constexpr bool operator != (const DynamicSet& other) const {
    return !(
      _items == other._items
    );
  }

  constexpr bool operator < (const DynamicSet& other) const {
    return _items < other._items;
  }

  constexpr bool operator > (const DynamicSet& other) const {
    return other._items < _items;
  }

  /*! 
   * Directly maps the container to a STL set without modifying the contained
   * elements
   */
  std::set<T> toSTL() const {
    std::set<T> returnSet;

    for(const auto& element : *this) {
      returnSet.insert(element);
    }

    return returnSet;
  }

  /*! 
   * Maps the contained elements to an STL set with a (possibly modifying) 
   * mapping function. This allows, e.g. a DynamicSet of DynamicArrays to be
   * directly mapped to an STL set of vectors.
   */
  template<typename MapFunction> 
  std::set<
    traits::functionReturnType<MapFunction, T>
  > mapToSTL(
    MapFunction&& function
  ) const {
    std::set<
      traits::functionReturnType<MapFunction, T>
    > returnSet;

    for(const auto& element : *this) {
      auto insertResultPair = returnSet.insert(
        function(element)
      );

      assert(insertResultPair.second);
    }

    return returnSet;
  }
};

//! Helper function to create a DynamicSet specifying only the maximum size
template<
  size_t nItems,
  typename T,
  template<typename, size_t> class ArrayType,
  class Comparator = std::less<T>
> constexpr DynamicSet<T, nItems, Comparator> makeDynamicSet(
  const ArrayType<T, nItems>& array
) {
  return DynamicSet<T, nItems, Comparator>(array);
}

} // namespace ConstexprMagic

#endif
