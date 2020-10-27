/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Compile-time container type traits
 *
 * Helper file that allows the use of container traits such as whether specific
 * member functions exist in SFINAE enable_ifs.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONTAINER_TRAITS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONTAINER_TRAITS_H

#include "Molassembler/Temple/Traits.h"

#include <utility>

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Traits {
namespace Detail {

/* All that is needed for boost's is_detected trickery (minus compatibility for
 * old compilers), replace as soon as minimum required boost version is >= 1.67
 */

template<class...>
struct make_void {
  using type = void;
};

template<class T>
using detector_t = typename make_void<T>::type;

template<class Default, class, template<class...> class, class...>
struct detector {
    using value_t = std::false_type;
    using type = Default;
};

template<class Default, template<class...> class Op, class... Args>
struct detector<Default, detector_t<Op<Args...> >, Op, Args...> {
    using value_t = std::true_type;
    using type = Op<Args...>;
};

struct nonesuch {};

template<template<class...> class Op, class... Args>
using is_detected = typename detector<nonesuch, void, Op, Args...>::value_t;

template<template<class...> class Op, class... Args>
constexpr bool is_detected_v = is_detected<Op, Args...>::value;

template<class Container>
using TestHasInsert = decltype(
  std::declval<Container>().insert(
    std::declval<getValueType<Container>>()
  )
);

template<class Container>
using TestHasPushBack = decltype(
  std::declval<Container>().push_back(
    std::declval<getValueType<Container>>()
  )
);

template<class Container>
using TestHasEmplace = decltype(
  std::declval<Container>().emplace(
    std::declval<getValueType<Container>>()
  )
);

template<class Container>
using TestHasEmplaceBack = decltype(
  std::declval<Container>().emplace_back(
    std::declval<getValueType<Container>>()
  )
);

template<class Container>
using TestHasSize = decltype(std::declval<Container>().size());

template<class Container>
using TestHasReserve = decltype(std::declval<Container>().reserve(0));

template<class Container>
using TestIsTuplelike = decltype(std::tuple_size<Container>::value);

template<class Container>
using TestIsPairlike = decltype(std::declval<Container>().first);

} // namespace Detail

template<class Container>
struct hasInsert : std::integral_constant<bool, Detail::is_detected_v<Detail::TestHasInsert, Container>> {};

template<class Container>
struct hasPushBack : std::integral_constant<bool, Detail::is_detected_v<Detail::TestHasPushBack, Container>> {};

template<class Container>
struct hasEmplace : std::integral_constant<bool, Detail::is_detected_v<Detail::TestHasEmplace, Container>> {};

template<class Container>
struct hasEmplaceBack : std::integral_constant<bool, Detail::is_detected_v<Detail::TestHasEmplaceBack, Container>> {};

template<class Container>
struct hasSize : std::integral_constant<bool, Detail::is_detected_v<Detail::TestHasSize, Container>> {};

template<class Container>
struct hasReserve : std::integral_constant<bool, Detail::is_detected_v<Detail::TestHasReserve, Container>> {};

template<class Container>
struct isTuplelike : std::integral_constant<bool, Detail::is_detected_v<Detail::TestIsTuplelike, Container>> {};

template<class Container>
struct isPairlike : std::integral_constant<bool, Detail::is_detected_v<Detail::TestIsPairlike, Container> && !isTuplelike<Container>::value> {};

} // namespace Traits
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
