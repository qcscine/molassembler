#include <type_traits>
#include <limits>
#include <exception>

#include "constexpr/Math.h"

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

  explicit Bitmask(EnumType a) {
    static_assert(
      std::is_unsigned<Underlying>::value
      && std::is_integral<Underlying>::value,
      "Underlying type for this enum type must unsigned and integral"
    );

    value = 1 << static_cast<Underlying>(a);
  }

  explicit Bitmask(Underlying a) {
    value = a;
  }

  Bitmask operator | (const EnumType& a) const {
    Underlying v = static_cast<Underlying>(a);

    if(v > maximum) {
      throw std::domain_error(
        "This enum has too many options to be representable as a bitmask "
        "in the specified underlying type."
      );
    }

    return Bitmask {
      value | (1 << v)
    };
  }

  inline bool isSet(const EnumType& a) const {
    return (
      value & (
        1 << static_cast<Underlying>(a)
      )
    ) > 0;
  }

  inline bool operator & (const EnumType& a) const {
    return isSet(a);
  }

  inline bool operator [] (const EnumType& a) const {
    return isSet(a);
  }
};

template<typename EnumType>
Bitmask<EnumType> make_bitmask(EnumType a) {
  return Bitmask<EnumType> {a};
}

} // namespace temple
