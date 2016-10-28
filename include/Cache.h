#include <experimental/optional>
#include <experimental/any>
#include <map>
#include <string>
#include <functional>
#include <cassert>
#include <vector>

/* TODO
 */

class Cache {
private:
/* Private members */
  std::map<
    std::string,
    std::experimental::any
  > _cache;
  std::map<
    std::string,
    std::function<
      std::experimental::any()
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
          std::experimental::any()  
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
  T getOrGenerate(const std::string& key) {
    if(_cache.count(key) == 1) {
      return std::experimental::any_cast<T>(
        _cache.at(key)
      );
    } else {
      // if this is false, user has violated contract
      assert(_generationMap.count(key) == 1); 
      _cache.emplace(
        key,
        _generationMap.at(key)() // calling it!
      );
      return std::experimental::any_cast<T>(
        _cache.at(key)
      );
    } 
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
  std::experimental::optional<T> getOption(const std::string& key) const {
    if(_cache.count(key) == 1) {
      return std::experimental::optional<T>(
        std::experimental::any_cast<T>(
          _cache.at(key)
        )
      );
    } else return std::experimental::nullopt;
  }


  bool hasValid(const std::string& key) const {
    return (_cache.count(key) > 0)
      ? true
      : false;
  }
};
