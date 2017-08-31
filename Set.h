#ifndef INCLUDE_CONSTEXPR_MAGIC_SET_H
#define INCLUDE_CONSTEXPR_MAGIC_SET_H

#include <array>

#include "Containers.h"
#include "Array.h"

/*! @file
 *
 * Minimal compile-time fixed-size set. Inserting elements returns a different
 * type signature.
 */

namespace ConstexprMagic {

namespace detail {

template<typename T>
struct defaultComparator {
  constexpr bool operator() (const T& a, const T& b) {
    return a < b;
  }
};

} // namespace detail

template<
  typename T,
  class Comparator = std::less<T>,
  size_t nItems = 0
> class Set {
private:
  /* C++17 -> switch to std::array, reference operator[] and reference at()
   * become constexpr, remember to change constIterator type name
   */
  using ArrayType = ConstexprMagic::Array<T, nItems>;

  ArrayType _items;
  Comparator _comparator;

public:
  constexpr Set() {}
  constexpr Set(const ArrayType& items) : _items(items) {}
  constexpr Set(ArrayType&& items) : _items(items) {}

  constexpr bool contains(const T& item) const {
    // C++17 turn into if constexpr-else
    if(nItems == 0) {
      return false;
    }

    auto bound = lowerBound<T, Comparator>(_items, item, _comparator);
    return !(item < _items.at(bound));
  }

  //! Insertion can involve up to 3N T moves.
  // TODO does not use Comparator, could be better -> see DynamicSet
  constexpr auto insert(const T& item) const {
    return Set<T, Comparator, nItems + 1>(
      insertIntoSorted(_items, item, _comparator)
    );
  }

  constexpr typename ArrayType::constIterator begin() const {
    return _items.begin();
  }

  constexpr typename ArrayType::constIterator end() const {
    return _items.end();
  }

  constexpr size_t size() const {
    return nItems;
  }

  constexpr bool operator == (const Set& other) const {
    return arraysEqual(_items, other._items);
  }

  constexpr bool operator != (const Set& other) const {
    return !arraysEqual(_items, other._items);
  }

  constexpr bool operator < (const Set& other) const {
    return arraysLess(_items, other._items);
  }

  constexpr bool operator > (const Set& other) const {
    return arraysLess(other._items, _items);
  }
};

template<
  typename T,
  class Comparator = std::less<T>,
  template<typename, size_t> class ArrayType,
  size_t nItems
> constexpr Set<T, Comparator, nItems> makeSetFromSortedArray(
  const ArrayType<T, nItems>& array
) {
  return Set<T, Comparator, nItems>(array);
}

} // namespace ConstexprMagic

#endif
