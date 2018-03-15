#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_MAP_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_MAP_H

#include "DynamicSet.h"
#include "Pair.h"

/*! @file
 *
 * Implements a constexpr container much like std::map, except with fixed
 * maximal size.
 */

namespace constable {

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

  struct OnlyFirstEquality {
    constexpr static auto _keyComparator = std::equal_to<KeyType>();

    constexpr OnlyFirstEquality() {};

    constexpr bool operator() (
      const PairType& a,
      const PairType& b
    ) const {
      return _keyComparator(a.first, b.first);
    }
  };

  using SetType = DynamicSet<PairType, N, OnlyFirstComparator, OnlyFirstEquality>;

  SetType _items;

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

  constexpr MappedType at(const KeyType& key) const {
    PairType pair {key, MappedType {}};
    auto keyOptional = _items.getOption(pair);
    
    if(!keyOptional.hasValue()) {
      throw "No such key in this DynamicMap!";
    }

    return keyOptional.value().second;
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

  using constIterator = typename SetType::constIterator;
  using const_iterator = constIterator;

  constexpr constIterator begin() const {
    return _items.begin();
  }

  constexpr constIterator end() const {
    return _items.end();
  }

  constexpr bool operator == (const DynamicMap& other) const {
    return _items == other._items;
  }

  constexpr bool operator != (const DynamicMap& other) const {
    return !(
      *this == other
    );
  }

  constexpr bool operator < (const DynamicMap& other) const {
    return _items < other._items;
  }

  constexpr bool operator > (const DynamicMap& other) const {
    return other._items < _items;
  }
};


} // namespace constable

#endif
