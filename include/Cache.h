#include <boost/optional.hpp>
#include <boost/any.hpp>
#include <map>
#include <string>
#include <functional>
#include <cassert>
#include <vector>

class Cache {
private:
/* Private members */
  std::map<
    std::string,
    boost::any
  > _cache;

  std::map<
    std::string,
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
        std::string,
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
  template<typename T>
  void add(const std::string& key, const T& value) {
    _cache.emplace(
      key,
      value // perfect forwarding
    );
  }


  template<typename T>
  T getGeneratable(const std::string& key) {
    // if this is false, user has violated contract
    assert(_generationMap.count(key) == 1);

    if(_cache.count(key) == 1) {
      return boost::any_cast<T>(
        _cache.at(key)
      );
    } else {
      _cache.emplace(
        key,
        _generationMap.at(key)() // calling it!
      );
      return boost::any_cast<T>(
        _cache.at(key)
      );
    } 
  }

  /* C++17
   * get rid of the raw pointer using map's extract and insert
   */
  template<typename T>
  void changeGeneratable(
    const std::string& key,
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

  void invalidate(const std::string& key) {
    if(_cache.count(key) == 1) {
      _cache.erase(key);
    }
  }

  /* Information */
  template<typename T>
  boost::optional<T> getOption(const std::string& key) const {
    if(_cache.count(key) == 1) {
      return boost::optional<T>(
        boost::any_cast<T>(
          _cache.at(key)
        )
      );
    } else return {};
  }


  bool has(const std::string& key) const {
    return (_cache.count(key) > 0)
      ? true
      : false;
  }
};
