#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_BITMASK_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_BITMASK_H

#include "temple/constexpr/Math.h"

#include <stdexcept>

/*!@file
 *
 * Contains a bitmask implementation for strong enums with explicit underlying
 * types and unmodified representational values.
 */

namespace temple {

template<typename EnumType>
struct Bitmask {
//!@name Types
//!@{
  using Underlying = std::underlying_type_t<EnumType>;
//!@}

  static constexpr Underlying maximum = temple::Math::floor(
    temple::Math::log(
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

  constexpr bool operator & (const EnumType& a) const {
    return isSet(a);
  }

  constexpr bool operator [] (const EnumType& a) const {
    return isSet(a);
  }
//!@}
};

template<typename EnumType>
constexpr Bitmask<EnumType> make_bitmask(EnumType a) {
  return Bitmask<EnumType> {a};
}

} // namespace temple

#endif
