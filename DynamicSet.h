#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_SET_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_SET_H

#include <set>

#include "DynamicArray.h"
#include "BTree.h"

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
  class LessThanPredicate = std::less<T>,
  class EqualityPredicate = std::equal_to<T>
> class DynamicSet {
private:
  using TreeType = BTree<T, 3, nItems, LessThanPredicate, EqualityPredicate>;

  TreeType _tree;

public:
  constexpr DynamicSet() {}

  //! Warning: These constructors expect ordered arrays!
  template<
    template<typename, size_t> class ArrayType,
    size_t size
  > constexpr DynamicSet(const ArrayType<T, size>& items) {
    for(const auto& item : items) {
      _tree.insert(item);
    }
  }

  //! Checks if an element exists. O(log N)
  constexpr bool contains(const T& item) const {
    return _tree.contains(item);
  }

  //! Insertion an element into the set. O(N)
  constexpr void insert(const T& item) {
    _tree.insert(item);
  }

  constexpr Optional<T> getOption(const T& item) const {
    return _tree.getOption(item);
  }

  constexpr void clear() {
    _tree.clear();
  }

  constexpr typename TreeType::constIterator begin() const {
    return _tree.begin();
  }

  constexpr typename TreeType::constIterator end() const {
    return _tree.end();
  }

  constexpr size_t size() const {
    return _tree.size();
  }

  constexpr bool operator == (const DynamicSet& other) const {
    return _tree == other._tree;
  }

  constexpr bool operator != (const DynamicSet& other) const {
    return !(
      _tree == other._tree
    );
  }

  constexpr bool operator < (const DynamicSet& other) const {
    return _tree < other._tree;
  }

  constexpr bool operator > (const DynamicSet& other) const {
    return other._tree < _tree;
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
  class LessThanPredicate = std::less<T>,
  class EqualityPredicate = std::equal_to<T>
> constexpr DynamicSet<T, nItems, LessThanPredicate> makeDynamicSet(
  const ArrayType<T, nItems>& array
) {
  return DynamicSet<T, nItems, LessThanPredicate, EqualityPredicate>(array);
}

} // namespace ConstexprMagic

#endif
