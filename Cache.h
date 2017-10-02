#ifndef INCLUDE_TEMPLATE_MAGIC_CACHE_H
#define INCLUDE_TEMPLATE_MAGIC_CACHE_H

#include <boost/any.hpp>
#include <boost/optional.hpp>

#include <map>
#include <string>
#include <functional>
#include <cassert>
#include <vector>

/*! @file
 *
 * Contains the Cache class implementation.
 */

namespace TemplateMagic {

/*! A minimal cache-like class.
 */
template<typename KeyType, typename ValueType>
class MinimalCache {
private:
/* Private members */
  //! Cache data
  std::map<
    KeyType,
    ValueType
  > _cache;

public:
/* Public member functions */
  /* Modification */
  //! Adds a data value for a key value into the cache.
  void add(const KeyType& key, const ValueType& value) {
    _cache.emplace(
      key,
      value
    );
  }

  //! Invalidates the entire cache, removing all stored data
  void invalidate() {
    _cache.clear();
  }

  //! Selectively invalidates cache entries
  void invalidate(const KeyType& key) {
    if(_cache.count(key) == 1) {
      _cache.erase(key);
    }
  }

  /* Information */
  //! Fetches a cache entry via an optional
  boost::optional<ValueType&> getOption(const KeyType& key) const {
    if(_cache.count(key) == 1) {
      return boost::optional<ValueType>(
        _cache.at(key)
      );
    } else return {};
  }

  //! Tests whether the cache contains an entry for a key
  bool has(const KeyType& key) const {
    return (_cache.count(key) > 0)
      ? true
      : false;
  }

  const ValueType& get(const KeyType& key) const {
    if(_cache.count(key) == 0) {
      throw "Fetching member in Cache whose key does not exist!";
    }

    return _cache.at(key);
  }
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
private:
/* Private members */
  //! Cache data
  std::map<
    KeyType,
    boost::any
  > _cache;

  //! Map of key values to generating functions
  std::map<
    KeyType,
    std::function<
      boost::any()
    >
  > _generationMap;

public:
/* Constructors */
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
      _generationMap[pair.first] = pair.second;
    }
  }


/* Public member functions */
  /* Modification */
  //! Adds a data value for a key value into the cache.
  template<typename T>
  void add(const KeyType& key, const T& value) {
    _cache.emplace(
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
    assert(_generationMap.count(key) == 1);

    if(_cache.count(key) == 1) {
      return boost::any_cast<T>(
        _cache.at(key)
      );
    } 

    _cache.emplace(
      key,
      _generationMap.at(key)() // calling it!
    );

    return boost::any_cast<T>(
      _cache.at(key)
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
    assert(_generationMap.count(key) == 1);
    
    // if the generatable does not exist yet, generate it
    if(_cache.count(key) == 0) {
      _cache.emplace(
        key,
        _generationMap.at(key)()
      );
    }

    // modify it
    modifyingUnaryFunction(
      boost::any_cast<T>(
        &(
          _cache.find(key) -> second
        ) 
      )
    );

  }

  //! Invalidates the entire cache, removing all stored data
  void invalidate() {
    _cache.clear();
  }

  //! Selectively invalidates cache entries
  void invalidate(const KeyType& key) {
    if(_cache.count(key) == 1) {
      _cache.erase(key);
    }
  }

  /* Information */
  //! Fetches a cache entry via an optional
  template<typename T>
  boost::optional<T> getOption(const KeyType& key) const {
    if(_cache.count(key) == 1) {
      return boost::optional<T>(
        boost::any_cast<T>(
          _cache.at(key)
        )
      );
    } else return {};
  }

  //! Tests whether the cache contains an entry for a key
  bool has(const KeyType& key) const {
    return (_cache.count(key) > 0)
      ? true
      : false;
  }
};

} //  namespace TemplateMagic

#endif
