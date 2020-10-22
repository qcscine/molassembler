/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Recursive serialization for debugging involving common STL containers.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_STRINGIFY_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_STRINGIFY_H

#include "Molassembler/Temple/Traits.h"

#include <set>
#include <unordered_set>
#include <list>
#include <forward_list>
#include <map>
#include <unordered_map>
#include <tuple>

#include "boost/optional.hpp"

namespace Scine {
namespace Molassembler {
namespace Temple {

/*!
 * Condenses an iterable container into a comma-separated string of string
 * representations of its contents. Requires container iterators to satisfy
 * ForwardIterators and the contained type to be a valid template
 * argument for std::to_string, which in essence means this works only for (a
 * few) STL containers and (most) built-in datatypes.
 */
template<class Container>
std::enable_if_t<
  !std::is_same<
    Traits::getValueType<Container>,
    std::string
  >::value,
  std::string
> condense(
  const Container& container,
  const std::string& joiningChar = ","
) {
  using namespace std::string_literals;

  std::string representation;

  for(auto it = std::begin(container); it != std::end(container); /*-*/) {
    representation += std::to_string(*it);
    if(++it != std::end(container)) {
      representation += joiningChar;
    }
  }

  return representation;
}

template<class Container> std::enable_if_t<
  std::is_same<
    Traits::getValueType<Container>,
    std::string
  >::value,
  std::string
> condense(
  const Container& container,
  const std::string& joiningChar = ","
) {
  using namespace std::string_literals;

  std::string representation;

  for(auto it = std::begin(container); it != std::end(container); /*-*/) {
    representation += *it;
    if(++it != std::end(container)) {
      representation += joiningChar;
    }
  }

  return representation;
}

namespace Detail {

template<class T, class R = void>
struct enable_if_type { using type = R; };

template<class T, class Enable = void>
struct isMapType : std::false_type {};

template<class T>
struct isMapType<T, typename enable_if_type<typename T::mapped_type>::type> : std::true_type {};

} // namespace Detail

using namespace std::string_literals;

inline std::string stringify(const std::string& x) {
  return x;
}

template<typename T>
std::string stringify(const boost::optional<T>& optional);

template<typename T>
std::enable_if_t<
  std::is_enum<T>::value,
  std::string
> stringify(const T& enumType);

template<typename T>
std::enable_if_t<
  std::is_arithmetic<T>::value,
  std::string
> stringify(const T& a);

template<typename T, typename U>
std::string stringify(const std::pair<T, U>& pair);

template<typename T>
std::string stringify(const std::vector<T>& vec);

template<typename T, std::size_t size>
std::string stringify(const std::array<T, size>& arr);

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

template<typename ... Ts>
std::string stringify(const std::tuple<Ts...>& tuple);

template<typename T>
std::enable_if_t<
  Detail::isMapType<T>::value,
  std::string
> stringifyMap(const T& map);

template<typename T>
std::string stringifyContainer(const T& container);

template<typename Container, class ElementStringifier>
std::string stringifyContainer(const Container& container, ElementStringifier&& stringifier);

template<typename TupleType, size_t ... Inds>
std::string stringifyTuple(
  const TupleType& tuple,
  std::index_sequence<Inds...> /* inds */
);


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
  return "("s + stringify(pair.first) + ", "s + stringify(pair.second) + ")"s;
}

template<typename T>
std::string stringify(const std::vector<T>& vec) {
  return stringifyContainer(vec);
}

template<typename T, size_t size>
std::string stringify(const std::array<T, size>& arr) {
  return stringifyContainer(arr);
}

template<typename T>
std::string stringify(const std::set<T>& set) {
  return stringifyContainer(set);
}

template<typename T>
std::string stringify(const std::unordered_set<T>& set) {
  return stringifyContainer(set);
}

template<typename T>
std::string stringify(const std::list<T>& list) {
  return stringifyContainer(list);
}

template<typename T, typename U>
std::string stringify(const std::map<T, U>& map) {
  return stringifyMap(map);
}

template<typename T, typename U>
std::string stringify(const std::unordered_map<T, U>& map) {
  return stringifyMap(map);
}

template<typename ... Ts>
std::string stringify(const std::tuple<Ts...>& tuple) {
  return stringifyTuple(tuple, std::make_index_sequence<sizeof...(Ts)> {});
}

template<typename T>
std::string stringify(const boost::optional<T>& optional) {
  if(optional) {
    return "Some "s + stringify(optional.value());
  }

  return "None";
}

template<typename T>
std::enable_if_t<
  std::is_enum<T>::value,
  std::string
> stringify(const T& enumType) {
  return "Enum="s + stringify(
    static_cast<
      std::underlying_type_t<T>
    >(enumType)
  );
}

template<typename T>
std::enable_if_t<
  Detail::isMapType<T>::value,
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
  std::string representation = "[";

  for(auto it = container.begin(); it != container.end(); /*-*/) {
    representation += stringify(*it);
    if(++it != container.end()) {
      representation += ", ";
    }
  }

  return representation + "]"s;
}

template<typename Container, class ElementStringifier>
std::string stringifyContainer(const Container& container, ElementStringifier&& stringifier) {
  std::string representation = "{";

  for(auto it = container.begin(); it != container.end(); /*-*/) {
    representation += stringifier(*it);
    if(++it != container.end()) {
      representation += ", ";
    }
  }

  return representation + "}"s;
}

template<typename TupleType, size_t ... Inds>
std::string stringifyTuple(
  const TupleType& tuple,
  std::index_sequence<Inds...> /* inds */
) {
  std::array<std::string, sizeof...(Inds)> individualStringifys {
    stringify(std::get<Inds>(tuple))...
  };

  return stringifyContainer(individualStringifys);
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
