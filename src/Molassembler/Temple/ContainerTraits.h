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
#include "boost/type_traits/is_detected.hpp"

#include <utility>

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Traits {
namespace Detail {

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

} // namespace Detail

template<class Container>
struct hasInsert : std::integral_constant<bool, boost::is_detected_v<Detail::TestHasInsert, Container>> {};

template<class Container>
struct hasPushBack : std::integral_constant<bool, boost::is_detected_v<Detail::TestHasPushBack, Container>> {};

template<class Container>
struct hasEmplace : std::integral_constant<bool, boost::is_detected_v<Detail::TestHasEmplace, Container>> {};

template<class Container>
struct hasEmplaceBack : std::integral_constant<bool, boost::is_detected_v<Detail::TestHasEmplaceBack, Container>> {};

template<class Container>
struct hasSize : std::integral_constant<bool, boost::is_detected_v<Detail::TestHasSize, Container>> {};

template<class Container>
struct hasReserve : std::integral_constant<bool, boost::is_detected_v<Detail::TestHasReserve, Container>> {};

} // namespace Traits
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
