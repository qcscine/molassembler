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

/*! Template base class for SFINAE expression validity checking
 *
 * @note The integer-long substitution trick is explained in Tricks.rst
 */
template<class>
struct sfinae_true : std::true_type{};

namespace Detail {

/* Short explanation: If the expression EXPR in sfinae_true<EXPR> returns a
 * valid type, sfinae_true can be instantiated, and the default-interpretation
 * of 0 as an int succeeds. If this fails, the backup-interpretation of 0 as a
 * long is chosen and returns std::false_type.
 */
template<class Container>
static auto testHasInsert(int)
  -> sfinae_true<
    decltype(
      std::declval<Container>().insert(
        std::declval<
          getValueType<Container>
        >()
      )
    )
  >;

template<class Container>
static auto testHasInsert(long) -> std::false_type;
} // namespace Detail

template<class Container>
struct hasInsert : decltype(Detail::testHasInsert<Container>(0)){};

namespace Detail {
template<class Container>
static auto testHasPushBack(int)
  -> sfinae_true<
    decltype(
      std::declval<Container>().push_back(
        std::declval<
          getValueType<Container>
        >()
      )
    )
  >;

template<class Container>
static auto testHasPushBack(long) -> std::false_type;
} // namespace Detail

template<class Container>
struct hasPushBack : decltype(Detail::testHasPushBack<Container>(0)){};

namespace Detail {
template<class Container>
static auto testHasEmplace(int)
  -> sfinae_true<
    decltype(
      std::declval<Container>().emplace(
        std::declval<
          getValueType<Container>
        >()
      )
    )
  >;

template<class Container>
static auto testHasEmplace(long) -> std::false_type;
} // namespace Detail

template<class Container>
struct hasEmplace : decltype(Detail::testHasEmplace<Container>(0)){};

namespace Detail {
template<class Container>
static auto testHasEmplaceBack(int)
  -> sfinae_true<
    decltype(
      std::declval<Container>().emplace_back(
        std::declval<
          getValueType<Container>
        >()
      )
    )
  >;

template<class Container>
static auto testHasEmplaceBack(long) -> std::false_type;
} // namespace Detail

template<class Container>
struct hasEmplaceBack : decltype(Detail::testHasEmplaceBack<Container>(0)){};


namespace Detail {
template<class Container>
static auto testHasSize(int)
  -> sfinae_true<
    decltype(
      std::declval<Container>().size()
    )
  >;

template<class Container>
static auto testHasSize(long) -> std::false_type;
} // namespace Detail

template<class Container>
struct hasSize : decltype(Detail::testHasSize<Container>(0)){};

namespace Detail {
template<class Container>
static auto testHasReserve(int)
  -> sfinae_true<
    decltype(
      std::declval<Container>().reserve()
    )
  >;

template<class Container>
static auto testHasReserve(long) -> std::false_type;
} // namespace Detail

template<class Container>
struct hasReserve : decltype(Detail::testHasReserve<Container>(0)){};

} // namespace Traits
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
