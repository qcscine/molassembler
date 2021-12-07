/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Helper functions for dealing with optional types
 *
 * Implements optional-returning function composition syntactic sugar to avoid
 * repetitive patterns when dealing with lots of optionals.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_OPTIONALS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_OPTIONALS_H

#include <boost/optional.hpp>
#include "Molassembler/Temple/Traits.h"

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Optionals {

/*! @brief Monadic bind with function of signature T -> U
 *
 * @tparam UnaryFunction: Function of signature T -> U
 * @tparam OptionalType: Optional type being used
 *
 * @returns OptionalType\<U\>
 */
template<
  template<typename> class OptionalType,
  typename T,
  class UnaryFunction
>
auto map(const OptionalType<T>& optional, UnaryFunction&& function) {
  using U = decltype(function(std::declval<T>()));

  if(optional) {
    return OptionalType<U> {function(optional.value())};
  }

  return OptionalType<U> {};
}

/*! @brief Monadic bind with function of signature T -> Optional\<U\>
 *
 * @tparam UnaryFunction: Function of signature T -> Optional\<U\>
 * @tparam OptionalType: Passed optional type
 *
 * @returns OptionalType\<U\>
 */
template<
  template<typename> class OptionalType,
  typename T,
  class UnaryFunction
>
auto flatMap(const OptionalType<T>& optional, UnaryFunction&& function) {
  using OptionalU = decltype(function(std::declval<T>()));

  if(optional) {
    return OptionalU {function(optional.value())};
  }

  return OptionalU {};
}

template<class Map, typename Key>
auto mapFind(const Map& map, const Key& key) {
  using Value = decltype(std::declval<Map>().begin()->second);
  using OptionalValue = boost::optional<Value>;

  const auto iter = map.find(key);
  if(iter != std::end(map)) {
    return OptionalValue(iter->second);
  }

  return OptionalValue(boost::none);
}

} // namespace Optionals
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
