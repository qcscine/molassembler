#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_SET_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_SET_H

#include <set>

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
  class LessThanPredicate = std::less<T>
> class DynamicSet {
private:
  /* C++17 -> switch to std::array, reference operator[] and reference at()
   * become constexpr, remember to change constIterator type name
   */
  using ArrayType = ConstexprMagic::DynamicArray<T, nItems>;

  ArrayType _items;
  LessThanPredicate _predicate;

  constexpr void _moveElementsRightUntil(
    typename ArrayType::iterator position
  ) {
    // Add the last element in the array onto the end
    _items.push_back(
      _items.back()
    );

    // Move everything up to and including the lower bound back
    auto rightIter = _items.end();
    rightIter -= 2;
    auto leftIter = rightIter - 1;

    while(rightIter != position) {
      *rightIter = std::move(*leftIter);

      --leftIter;
      --rightIter;
    }
  }

public:
  constexpr DynamicSet() {}

  //! Warning: These constructors expect ordered arrays!
  constexpr DynamicSet(const ArrayType& items) : _items(items) {}
  constexpr DynamicSet(ArrayType&& items) : _items(items) {}

  //! Checks if an element exists. O(log N)
  constexpr bool contains(const T& item) const {
    if(_items.size() == 0) {
      return false;
    }

    return binarySearch(_items, item, _predicate) != _items.end();
  }

  /*!
   * Returns an iterator to the sought element, or the end iterator if the
   * element is not found. O(log N)
   */
  constexpr auto find(const T& item) const {
    return binarySearch(_items, item, _predicate);
  }

  /*!
   * Returns an iterator to the first element that is not smaller than the passed
   * value. O(log N)
   */
  constexpr typename DynamicArray<T, nItems>::iterator getLowerBound(const T& item) {
    return lowerBound<T, LessThanPredicate>(
      _items.begin(),
      _items.end(),
      item,
      _predicate
    );
  }

  /*!
   * Returns an iterator to the first element that is not smaller than the passed
   * value.
   *
   * O(log N)
   */
  constexpr typename DynamicArray<T, nItems>::const_iterator getLowerBound(const T& item) const {
    return lowerBound<T, LessThanPredicate>(
      _items.begin(),
      _items.end(),
      item,
      _predicate
    );
  }

  //! Inserts an element at a specified position. O(N)
  constexpr void insertAt(
    typename ArrayType::iterator insertPositionIter,
    const T& item
  ) {
    // In case there is no object larger than the one being inserted, add to end
    if(insertPositionIter == _items.end()) {
      _items.push_back(item);
    } else {
      _moveElementsRightUntil(insertPositionIter);

      // Copy in the item
      *insertPositionIter = item;
    }
  }

  //! Inserts an element at a specified position O(N)
  constexpr void insertAt(
    typename ArrayType::iterator insertPositionIter,
    T&& item
  ) {
    // In case there is no object larger than the one being inserted, add to end
    if(insertPositionIter == _items.end()) {
      _items.push_back(item);
    } else {
      _moveElementsRightUntil(insertPositionIter);

      // Move in the item
      *insertPositionIter = std::move(item);
    }
  }

  //! Insertion an element into the set. O(N)
  constexpr void insert(const T& item) {
    insertAt(
      getLowerBound(item),
      item
    );
  }

  //! Inserts an element into the set. O(N)
  constexpr void insert(T&& item) {
    insertAt(
      getLowerBound(item),
      std::forward<T>(item)
    );
  }

  constexpr bool lowerBoundMeansContains(
    const typename DynamicArray<T, nItems>::iterator& lowerBound,
    const T& item
  ) {
    if(_predicate(item, *lowerBound) || lowerBound == _items.end()) {
      return false;
    }

    return true;
  }

  constexpr bool lowerBoundMeansContains(
    const typename DynamicArray<T, nItems>::const_iterator& lowerBound,
    const T& item
  ) const {
    if(_predicate(item, *lowerBound) || lowerBound == _items.end()) {
      return false;
    }

    return true;
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
  class LessThanPredicate = std::less<T>
> constexpr DynamicSet<T, nItems, LessThanPredicate> makeDynamicSet(
  const ArrayType<T, nItems>& array
) {
  return DynamicSet<T, nItems, LessThanPredicate>(array);
}

} // namespace ConstexprMagic

#endif
