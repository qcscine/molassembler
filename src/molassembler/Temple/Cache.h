/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Cache classes
 *
 * Contains a map-like MinimalCache and a full-blown boost::any-Cache.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CACHE_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CACHE_H

#include <boost/any.hpp>
#include <boost/optional.hpp>

#include <map>
#include <string>
#include <functional>
#include <cassert>
#include <vector>

namespace Scine {
namespace temple {

//! A minimal wrapper around a map class with cache semantics
template<typename KeyType, typename ValueType>
class MinimalCache {
public:
  using key_type = KeyType;
  using value_type = ValueType;

//!@name Modification
//!@{
  //! Adds a data value for a key value into the cache.
  void add(const KeyType key, const ValueType value) {
    cache_.emplace(
      std::move(key),
      std::move(value)
    );
  }

  //! Invalidates the entire cache, removing all stored data
  void invalidate() {
    cache_.clear();
  }

  //! Selectively invalidates cache entries
  void invalidate(const KeyType& key) {
    if(cache_.count(key) == 1) {
      cache_.erase(key);
    }
  }
//!@}

//!@name Information
//!@{
  //! Fetches a cache entry via an optional
  boost::optional<const ValueType&> getOption(const KeyType& key) const {
    if(cache_.count(key) == 1) {
      return boost::optional<const ValueType&>(
        cache_.at(key)
      );
    }

    return {};
  }

  //! Tests whether the cache contains an entry for a key
  bool has(const KeyType& key) const {
    return cache_.count(key) > 0;
  }

  //! Fetches a particular key, assuming that it is currently cached
  const ValueType& get(const KeyType& key) const {
    if(cache_.count(key) == 0) {
      throw "Fetching member in Cache whose key does not exist!";
    }

    return cache_.at(key);
  }
//!@}

private:
/* Private members */
  //! Cache data
  std::map<KeyType, ValueType> cache_;
};

/*!
 *
 * A class that can store completely variant types in a type safe manner using
 * boost::any. Only the key type is fixed. Selective invalidation of cached data
 * is possible, and cache elements can also be set as generatable, i.e. a
 * function that generates the cache data can be passed. On first use, the cache
 * element is then generated.
 */
template<typename KeyType>
class Cache {
public:
//!@name Constructors
//!@{
  Cache() = default;
  Cache(
    const std::initializer_list<
      std::pair<
        KeyType,
        std::function<
          boost::any()
        >
      >
    >& initList
  ) {
    // add all generators
    for(const auto& pair: initList) {
      generationMap_[pair.first] = pair.second;
    }
  }
//!@}


//!@name Modification
//!@{
  //! Adds a data value for a key value into the cache.
  template<typename T>
  void add(const KeyType& key, const T& value) {
    cache_.emplace(
      key,
      value
    );
  }

  /*!
   * Fetches the value of a generatable cache item. Handles all corner cases of
   * whether the cache element is present or not.
   */
  template<typename T>
  T getGeneratable(const KeyType& key) {
    // if this is false, user has violated contract
    assert(generationMap_.count(key) == 1);

    if(cache_.count(key) == 1) {
      return boost::any_cast<T>(
        cache_.at(key)
      );
    }

    cache_.emplace(
      key,
      generationMap_.at(key)() // calling it!
    );

    return boost::any_cast<T>(
      cache_.at(key)
    );
  }

  /* C++17
   * get rid of the raw pointer using map's extract and insert
   */
  //! Permits the modification of a generatable cache item via a raw pointer.
  template<typename T>
  void changeGeneratable(
    const KeyType& key,
    std::function<
      void(T*)
    > modifyingUnaryFunction
  ) {
    assert(generationMap_.count(key) == 1);

    // if the generatable does not exist yet, generate it
    if(cache_.count(key) == 0) {
      cache_.emplace(
        key,
        generationMap_.at(key)()
      );
    }

    // modify it
    modifyingUnaryFunction(
      boost::any_cast<T>(
        &(
          cache_.find(key) -> second
        )
      )
    );

  }

  //! Invalidates the entire cache, removing all stored data
  void invalidate() {
    cache_.clear();
  }

  //! Selectively invalidates cache entries
  void invalidate(const KeyType& key) {
    if(cache_.count(key) == 1) {
      cache_.erase(key);
    }
  }
//!@}

//!@name Information
//!@{
  //! Fetches a cache entry via an optional
  template<typename T>
  boost::optional<T> getOption(const KeyType& key) const {
    if(cache_.count(key) == 1) {
      return boost::optional<T>(
        boost::any_cast<T>(
          cache_.at(key)
        )
      );
    }

    return {};
  }

  //! Tests whether the cache contains an entry for a key
  bool has(const KeyType& key) const {
    return (cache_.count(key) > 0);
  }
//!@}

private:
/* Private members */
  //! Cache data
  std::map<
    KeyType,
    boost::any
  > cache_;

  //! Map of key values to generating functions
  std::map<
    KeyType,
    std::function<
      boost::any()
    >
  > generationMap_;
};

} //  namespace temple
} // namespace Scine

#endif
