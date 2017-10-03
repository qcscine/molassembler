#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_SET_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_SET_H

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

  using constIterator = typename TreeType::constIterator;
  using const_iterator = constIterator;

  constexpr constIterator begin() const {
    return _tree.begin();
  }

  constexpr constIterator end() const {
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
