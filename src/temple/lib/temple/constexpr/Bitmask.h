#include <type_traits>
#include <limits>
#include <exception>

#include "Math.h"

namespace temple {

template<typename EnumType>
struct Bitmask {
  using Underlying = std::underlying_type_t<EnumType>;

  static constexpr Underlying maximum = temple::Math::floor(
    temple::Math::log(
      static_cast<double>(std::numeric_limits<Underlying>::max()),
      2.0
    )
  );

  Underlying value;

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

  constexpr Bitmask operator | (const EnumType& a) const {
    Underlying v = static_cast<Underlying>(a);

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

  inline constexpr bool isSet(const EnumType& a) const {
    return (
      value & (
        static_cast<Underlying>(1) << static_cast<Underlying>(a)
      )
    ) > 0;
  }

  inline constexpr bool operator & (const EnumType& a) const {
    return isSet(a);
  }

  inline constexpr bool operator [] (const EnumType& a) const {
    return isSet(a);
  }
};

template<typename EnumType>
constexpr Bitmask<EnumType> make_bitmask(EnumType a) {
  return Bitmask<EnumType> {a};
}

} // namespace temple
