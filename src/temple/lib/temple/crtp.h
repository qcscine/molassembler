#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CRTP_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CRTP_H

namespace temple {

namespace crtp {

//! Supplies the inequality operator from an implemented equality operator
template<typename T>
struct Inequality {
  constexpr bool operator != (const Inequality<T>& other) const {
    const auto& lhs = static_cast<const T&>(*this);
    const auto& rhs = static_cast<const T&>(other);

    return !(lhs == rhs);
  }
};

//! Supplies all operators from an implemented less-than operator
template <typename T>
struct AllOperators {
  static constexpr const T& getDerived(const AllOperators& base) {
    return static_cast<const T&>(base);
  }

  constexpr bool operator == (const AllOperators& other) const {
    auto lhs = getDerived(*this), rhs = getDerived(other);

    return !(lhs < rhs) && !(rhs < lhs);
  }

  constexpr bool operator != (const AllOperators& other) const {
    auto lhs = getDerived(*this), rhs = getDerived(other);

    return (lhs < rhs) || (rhs < lhs);
  }

  constexpr bool operator > (const AllOperators& other) const {
    auto lhs = getDerived(*this), rhs = getDerived(other);

    return rhs < lhs;
  }

  constexpr bool operator <= (const AllOperators& other) const {
    auto lhs = getDerived(*this), rhs = getDerived(other);

    return !(rhs < lhs);
  }

  constexpr bool operator >= (const AllOperators& other) const {
    auto lhs = getDerived(*this), rhs = getDerived(other);

    return !(lhs < rhs);
  }
};


} // namespace crtp

} // namespace temple

#endif
