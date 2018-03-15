#ifndef INCLUDE_TEMPLATE_MAGIC_STRINGIFY_H
#define INCLUDE_TEMPLATE_MAGIC_STRINGIFY_H

#include "Traits.h"
#include <set>
#include <unordered_set>
#include <list>
#include <forward_list>
#include <map>
#include <unordered_map>

/*! @file
 *
 * Recursive serialization for quicker debugging of common STL containers.
 */

namespace temple {

namespace detail {

template<class T, class R = void>  
struct enable_if_type { typedef R type; };

template<class T, class Enable = void>
struct isMapType : std::false_type {};

template<class T>
struct isMapType<T, typename enable_if_type<typename T::mapped_type>::type> : std::true_type {};

} // namespace detail

using namespace std::string_literals;

inline std::string stringify(const std::string& x) {
  return x;
}

template<typename T>
std::enable_if_t<
  std::is_arithmetic<T>::value,
  std::string
> stringify(const T& a);

template<typename T, typename U>
std::string stringify(const std::pair<T, U>& pair);

template<typename T>
std::string stringify(const std::vector<T>& vec);

template<typename T>
std::string stringify(const std::set<T>& set);

template<typename T>
std::string stringify(const std::unordered_set<T>& set);

template<typename T>
std::string stringify(const std::list<T>& list);

template<typename T>
std::string stringify(const std::forward_list<T>& list);

template<typename T, typename U>
std::string stringify(const std::map<T, U>& map);

template<typename T, typename U>
std::string stringify(const std::unordered_map<T, U>& map);

template<typename T>
std::enable_if_t<
  detail::isMapType<T>::value,
  std::string
> stringifyMap(const T& map);

template<typename T>
std::string stringifyContainer(const T& container);


// Definitions

template<typename T>
std::enable_if_t<
  std::is_arithmetic<T>::value,
  std::string
> stringify(const T& a) {
  return std::to_string(a);
}

template<typename T, typename U>
std::string stringify(const std::pair<T, U>& pair) {
  return "pair {"s + stringify(pair.first) + ", "s + stringify(pair.second) + "}"s;
}

template<typename T>
std::string stringify(const std::vector<T>& vec) {
  return "vector "s + stringifyContainer(vec);
}

template<typename T>
std::string stringify(const std::set<T>& set) {
  return "set "s + stringifyContainer(set);
}

template<typename T>
std::string stringify(const std::unordered_set<T>& set) {
  return "unord. set "s + stringifyContainer(set);
}

template<typename T>
std::string stringify(const std::list<T>& list) {
  return "list "s + stringifyContainer(list);
}

template<typename T>
std::string stringify(const std::forward_list<T>& list) {
  return "fwd. list "s + stringifyContainer(list);
}

template<typename T, typename U>
std::string stringify(const std::map<T, U>& map) {
  return "map "s + stringifyMap(map);
}

template<typename T, typename U>
std::string stringify(const std::unordered_map<T, U>& map) {
  return "unord. map "s + stringifyMap(map);
}

template<typename T>
std::enable_if_t<
  detail::isMapType<T>::value,
  std::string
> stringifyMap(const T& map) {
  std::string representation = "{";

  for(auto it = map.begin(); it != map.end(); /*-*/) {
    representation += stringify(it->first) + " -> "s + stringify(it->second);
    if(++it != map.end()) {
      representation += ", ";
    }
  }

  return representation + "}"s;
}

template<typename T>
std::string stringifyContainer(const T& container) {
  std::string representation = "{";

  for(auto it = container.begin(); it != container.end(); /*-*/) {
    representation += stringify(*it);
    if(++it != container.end()) {
      representation += ", ";
    }
  }

  return representation + "}"s;
}

} // namespace temple

#endif
