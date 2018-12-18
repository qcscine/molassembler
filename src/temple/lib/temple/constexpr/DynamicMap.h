/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief A std::map-like class based on BTree (but max size is space allocated)
 *
 * Implements a constexpr container much like std::map, except with fixed
 * maximal size.
 */

#ifndef INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_DYNAMIC_MAP_H
#define INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_DYNAMIC_MAP_H

#include "temple/constexpr/DynamicSet.h"
#include "temple/constexpr/Pair.h"

namespace temple {

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

    constexpr OnlyFirstComparator() = default;

    constexpr bool operator() (
      const PairType& a,
      const PairType& b
    ) const {
      return _keyComparator(a.first, b.first);
    }
  };

  struct OnlyFirstEquality {
    constexpr static auto _keyComparator = std::equal_to<KeyType>();

    constexpr OnlyFirstEquality() = default;

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
//!@name Special member functions
//!@{
  constexpr DynamicMap() = default;

  constexpr DynamicMap(DynamicMap&& other) noexcept
    : _items(other._items)
  {}

  constexpr DynamicMap(const DynamicMap& other)
    : _items(other._items)
  {}

  constexpr DynamicMap& operator = (const DynamicMap& other) {
    _items = other._items;
    return *this;
  }

  constexpr DynamicMap& operator = (DynamicMap&& other) noexcept {
    _items = other._items;
    return *this;
  }
//!@}

//!@name Modification
//!@{
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
      // searchIter is a const_iterator unfortunately, so need to go via index
      auto indexOfElement = static_cast<unsigned>(
        searchIter - _items.begin()
      );

      // Overwrite pair with same key
      *(_items.begin() + indexOfElement) = pair;
    }
  }

  constexpr void clear() {
    _items.clear();
  }
//!@}

//!@name Information
//!@{
  constexpr MappedType at(const KeyType& key) const {
    PairType pair {key, MappedType {}};
    auto keyOptional = _items.getOption(pair);

    if(!keyOptional.hasValue()) {
      throw "No such key in this DynamicMap!";
    }

    return keyOptional.value().second;
  }


  constexpr unsigned size() const {
    return _items.size();
  }
//!@}

//!@name Iterators
//!@{
  using const_iterator = typename SetType::const_iterator;

  constexpr const_iterator begin() const {
    return _items.begin();
  }

  constexpr const_iterator end() const {
    return _items.end();
  }
//!@}

//!@name Operators
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
//!@}
};


} // namespace temple

#endif
