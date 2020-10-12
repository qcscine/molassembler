/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief bitmask class
 *
 * Contains a bitmask implementation for strong enums with explicit underlying
 * types and unmodified representational values.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_BITMASK_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_BITMASK_H

#include "Molassembler/Temple/constexpr/Math.h"

#include <stdexcept>

namespace Scine {
namespace Molassembler {
namespace Temple {

/**
 * @brief
 *
 * @tparam EnumType Enum on which the bitmask should act. Requires the enum to
 * have strictly incrementing unsigned representation.
 */
template<typename EnumType>
struct Bitmask {
//!@name Types
//!@{
  using Underlying = std::underlying_type_t<EnumType>;
//!@}

  static constexpr Underlying maximum = Temple::Math::floor(
    Temple::Math::log(
      static_cast<double>(std::numeric_limits<Underlying>::max()),
      2.0
    )
  );

//!@name Public state
//!@{
  Underlying value;
//!@}

//!@name Constructors
//!@{
  explicit constexpr Bitmask() : value {0} {}

  explicit constexpr Bitmask(EnumType a) : value {
    static_cast<Underlying>(1) << static_cast<Underlying>(a)
  } {
    static_assert(
      std::is_unsigned<Underlying>::value
      && std::is_integral<Underlying>::value,
      "Underlying type for this enum type must unsigned and integral"
    );
  }

  explicit constexpr Bitmask(Underlying a) : value {a} {}
//!@}

//!@name Information
//!@{
  /*! @brief Checks whether an enum value is set in the bitmask
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr bool isSet(const EnumType& a) const {
    return (
      value & (
        static_cast<Underlying>(1) << static_cast<Underlying>(a)
      )
    ) > 0;
  }
//!@}

//!@name Operators
//!@{
  /*! @brief Create a new bitmask that also sets a particular enum value
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr Bitmask operator | (const EnumType& a) const {
    auto v = static_cast<Underlying>(a);

    if(v > maximum) {
      throw std::domain_error(
        "This enum has too many options to be representable as a bitmask "
        "in the specified underlying type."
      );
    }

    return Bitmask {
      value | (static_cast<Underlying>(1) << v)
    };
  }

  /*! @brief Set a particular enum value in this bitmask
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr void operator |= (const EnumType& a) {
    auto v = static_cast<Underlying>(a);

    if(v > maximum) {
      throw std::domain_error(
        "This enum has too many options to be representable as a bitmask "
        "in the specified underlying type."
      );
    }

    value = value | (static_cast<Underlying>(1) << v);
  }

  //! Check whether a particular enum value is set
  constexpr bool operator & (const EnumType& a) const {
    return isSet(a);
  }

  //! Check whether a particular enum value is set
  constexpr bool operator [] (const EnumType& a) const {
    return isSet(a);
  }
//!@}
};

template<typename EnumType>
constexpr Bitmask<EnumType> make_bitmask(EnumType a) {
  return Bitmask<EnumType> {a};
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
