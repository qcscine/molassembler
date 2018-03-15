#ifndef INCLUDE_MEMBER_FUNCTION_CACHE_H
#define INCLUDE_MEMBER_FUNCTION_CACHE_H

#include <boost/optional.hpp>
#include <boost/any.hpp>
#include <map>
#include <functional>
#include <cassert>

/*! @file
 *
 * Contains the implementation of a cache that facilitates the caching of the
 * result of class instance member functions.
 */

/*!
 * Facilitates the caching of the result of class instance member functions.
 */
template<typename KeyType, typename BaseType>
class MemFnCache {
public:
  //! The type signature of accepted functions for the generation of cache data
  using LambdaType  = std::function<
    boost::any(const BaseType&)
  >;

private:
/* Private members */
  std::map<KeyType, boost::any> _cache;
  std::map<KeyType, LambdaType> _generationMap;

public:
/* Constructors */
  MemFnCache() = delete;
  MemFnCache(
    const std::initializer_list<
      std::pair<KeyType, LambdaType>
    >& initList
  ) {
    // add all generators
    for(const auto& pair: initList) {
      _generationMap[pair.first] = pair.second;
    }
  }


/* Public member functions */
  /* Modification */
  template<typename T>
  void add(const KeyType& key, const T& value) {
    _cache.emplace(
      key,
      value // perfect forwarding
    );
  }


  template<typename T>
  T getGeneratable(const KeyType& key, const BaseType& baseRef) {
    // if this is false, user has violated contract
    assert(_generationMap.count(key) == 1);

    if(_cache.count(key) == 1) {
      return boost::any_cast<T>(
        _cache.at(key)
      );
    } 

    _cache.emplace(
      key,
      _generationMap.at(key)(baseRef) // calling it!
    );

    return boost::any_cast<T>(
      _cache.at(key)
    );
  }

  /* C++17
   * get rid of the raw pointer using map's extract and insert
   */
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

  void invalidate() {
    _cache.clear();
  }

  void invalidate(const KeyType& key) {
    if(_cache.count(key) == 1) {
      _cache.erase(key);
    }
  }

  /* Information */
  template<typename T>
  boost::optional<T> getOption(const KeyType& key) const {
    if(_cache.count(key) == 1) {
      return boost::optional<T>(
        boost::any_cast<T>(
          _cache.at(key)
        )
      );
    } else return boost::none;
  }


  bool has(const KeyType& key) const {
    if(_cache.count(key) > 0) {
      return true;
    }

    return false;
  }
};

#endif
