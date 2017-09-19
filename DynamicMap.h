#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_MAP_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_MAP_H

#include "DynamicSet.h"
#include "Pair.h"

/*! @file
 *
 * Implements a constexpr container much like std::map, except with fixed
 * maximal size.
 */

namespace ConstexprMagic {

/*!
 * A constexpr associative container with reasonably fast key-based lookup.
 *
 * Requires that MappedType is default-constructible
 */
template<typename KeyType, typename MappedType, size_t N>
class DynamicMap {
private:
  /* std::pair is NOT constexpr move or copy assignable, and this does not
   * appear to change with C++17
   */
  using PairType = Pair<KeyType, MappedType>;

  struct OnlyFirstComparator {
    constexpr static auto _keyComparator = std::less<KeyType>();

    constexpr OnlyFirstComparator() {};

    constexpr bool operator() (
      const PairType& a,
      const PairType& b
    ) const {
      return _keyComparator(a.first, b.first);
    }
  };

  DynamicSet<PairType, N, OnlyFirstComparator> _items;

public:
  constexpr DynamicMap() {}

  constexpr DynamicMap(DynamicMap&& other)
    : _items(other._items)
  {}

  constexpr DynamicMap(const DynamicMap& other)
    : _items(other._items)
  {}

  constexpr DynamicMap& operator = (const DynamicMap& other) {
    _items = other._items;
    return *this;
  }

  constexpr DynamicMap& operator = (DynamicMap&& other) {
    _items = other._items;
    return *this;
  }

  constexpr const MappedType& at(const KeyType& key) const {
    auto searchIter = _items.find(PairType {key, MappedType {}});

    if(searchIter == _items.end()) {
      throw "No such key in this DynamicMap!";
    }

    else return (*searchIter).second;
  }
      
  constexpr void insert(
    KeyType key,
    MappedType item
  ) {
    PairType pair { std::move(key), std::move(item) };

    if(_items.contains(pair)) {
      throw "Map already contains an item for this key!";
    }

    _items.insert(pair);
  }

  constexpr void insertOrUpdate(
    KeyType key,
    MappedType item
  ) {
    PairType pair { std::move(key), std::move(item) };

    auto searchIter = _items.find(pair);

    if(searchIter == _items.end()) {
      _items.insert(pair);
    } else {
      // searchIter is a constIterator unfortunately, so need to go via index
      unsigned indexOfElement = static_cast<unsigned>(
        searchIter - _items.begin()
      );

      // Overwrite pair with same key
      *(_items.begin() + indexOfElement) = pair;
    }
  }

  constexpr void clear() {
    _items.clear();
  }

  constexpr unsigned size() const {
    return _items.size();
  }
};


} // namespace ConstexprMagic

#endif
