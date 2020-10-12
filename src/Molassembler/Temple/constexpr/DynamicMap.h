/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief A std::map-like class based on BTree (but max size is space allocated)
 *
 * Implements a constexpr container much like std::map, except with fixed
 * maximal size.
 */

#ifndef INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_DYNAMIC_MAP_H
#define INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_DYNAMIC_MAP_H

#include "Molassembler/Temple/constexpr/DynamicSet.h"
#include "Molassembler/Temple/constexpr/Pair.h"

namespace Scine {
namespace Molassembler {
namespace Temple {

/*! @brief A constexpr associative container with reasonably fast key-based lookup.
 *
 * Requires that MappedType is default-constructible
 *
 * @tparam KeyType Type of the key in the map
 * @tparam MappedType Value type of the map
 * @tparam N Maximum number of elements in the map
 */
template<typename KeyType, typename MappedType, size_t N>
class DynamicMap {
private:
  /* std::pair is NOT constexpr move or copy assignable, and this does not
   * appear to change with C++17
   */
  using PairType = Pair<KeyType, MappedType>;

  /* A map is just a tree set where each node stores both a key and a value.
   * Comparisons and equality merely reference the keys. Value type lookup by
   * key is just searching for a node by key and returning the stored value
   * type there.
   */
  struct OnlyFirstComparator {
    constexpr static auto keyComparator_ = std::less<KeyType>();

    constexpr OnlyFirstComparator() = default;

    constexpr bool operator() (
      const PairType& a,
      const PairType& b
    ) const {
      return keyComparator_(a.first, b.first);
    }
  };

  struct OnlyFirstEquality {
    constexpr static auto keyComparator_ = std::equal_to<KeyType>();

    constexpr OnlyFirstEquality() = default;

    constexpr bool operator() (
      const PairType& a,
      const PairType& b
    ) const {
      return keyComparator_(a.first, b.first);
    }
  };

  using SetType = DynamicSet<PairType, N, OnlyFirstComparator, OnlyFirstEquality>;

  SetType items_;

public:
//!@name Special member functions
//!@{
  constexpr DynamicMap() = default;

  constexpr DynamicMap(DynamicMap&& other) noexcept
    : items_(other.items_)
  {}

  constexpr DynamicMap(const DynamicMap& other)
    : items_(other.items_)
  {}

  constexpr DynamicMap& operator = (const DynamicMap& other) {
    items_ = other.items_;
    return *this;
  }

  constexpr DynamicMap& operator = (DynamicMap&& other) noexcept {
    items_ = other.items_;
    return *this;
  }
//!@}

//!@name Modification
//!@{
  /*! @brief Inserts an element into the map
   *
   * @complexity{@math{\Theta(N \log N)}}
   */
  constexpr void insert(
    KeyType key,
    MappedType item
  ) {
    PairType pair { std::move(key), std::move(item) };

    if(items_.contains(pair)) {
      throw "Map already contains an item for this key!";
    }

    items_.insert(pair);
  }

  /*! @brief Inserts a key-value pair into the map or updates the mapped value
   *   if the key exists
   *
   * @complexity{@math{\Theta(N \log N)}}
   */
  constexpr void insertOrUpdate(
    KeyType key,
    MappedType item
  ) {
    PairType pair { std::move(key), std::move(item) };

    auto searchIter = items_.find(pair);

    if(searchIter == items_.end()) {
      items_.insert(pair);
    } else {
      // searchIter is a const_iterator unfortunately, so need to go via index
      auto indexOfElement = static_cast<unsigned>(
        searchIter - items_.begin()
      );

      // Overwrite pair with same key
      *(items_.begin() + indexOfElement) = pair;
    }
  }

  constexpr void clear() {
    items_.clear();
  }
//!@}

//!@name Information
//!@{
  /*! @brief Value lookup.
   *
   * @complexity{@math{\Theta(N \log N)}}
   */
  constexpr const MappedType& at(const KeyType& key) const {
    PairType pair {key, MappedType {}};
    auto keyOptional = items_.getOption(pair);

    if(!keyOptional.hasValue()) {
      throw "No such key in this DynamicMap!";
    }

    return keyOptional.value().second;
  }


  /*! @brief Number of items in the map
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr unsigned size() const {
    return items_.size();
  }
//!@}

//!@name Iterators
//!@{
  using const_iterator = typename SetType::const_iterator;

  constexpr const_iterator begin() const {
    return items_.begin();
  }

  constexpr const_iterator end() const {
    return items_.end();
  }
//!@}

//!@name Operators
//!@{
  constexpr bool operator == (const DynamicMap& other) const {
    return items_ == other.items_;
  }

  constexpr bool operator != (const DynamicMap& other) const {
    return !(*this == other);
  }

  constexpr bool operator < (const DynamicMap& other) const {
    return items_ < other.items_;
  }

  constexpr bool operator > (const DynamicMap& other) const {
    return other.items_ < items_;
  }
//!@}
};

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
