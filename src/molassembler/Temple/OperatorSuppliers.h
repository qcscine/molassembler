/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Operator-supplying CRTP base classes
 *
 * Permits the provision of more operators from fewer initial operator
 * implementations.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_OPERATOR_SUPPLIERS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_OPERATOR_SUPPLIERS_H

namespace Scine {
namespace temple {

//! @brief Curiously recurring template pattern classes
namespace crtp {

//! Supplies the inequality operator from an implemented equality operator
template<typename T>
struct InequalityFromEquality {
  constexpr bool operator != (const InequalityFromEquality<T>& other) const {
    const auto& lhs = static_cast<const T&>(*this);
    const auto& rhs = static_cast<const T&>(other);

    return !(lhs == rhs);
  }
};

//! Supplies all operators from an implemented less-than operator
template <typename T>
struct AllOperatorsFromLessThan {
  static constexpr const T& getDerived(const AllOperatorsFromLessThan& base) {
    return static_cast<const T&>(base);
  }

  constexpr bool operator == (const AllOperatorsFromLessThan& other) const {
    const T& lhs = getDerived(*this), rhs = getDerived(other);

    return !(lhs < rhs) && !(rhs < lhs);
  }

  constexpr bool operator != (const AllOperatorsFromLessThan& other) const {
    const T& lhs = getDerived(*this), rhs = getDerived(other);

    return (lhs < rhs) || (rhs < lhs);
  }

  constexpr bool operator > (const AllOperatorsFromLessThan& other) const {
    const T& lhs = getDerived(*this), rhs = getDerived(other);

    return rhs < lhs;
  }

  constexpr bool operator <= (const AllOperatorsFromLessThan& other) const {
    const T& lhs = getDerived(*this), rhs = getDerived(other);

    return !(rhs < lhs);
  }

  constexpr bool operator >= (const AllOperatorsFromLessThan& other) const {
    const T& lhs = getDerived(*this), rhs = getDerived(other);

    return !(lhs < rhs);
  }
};

/*! @brief Generates all operators using a method returning a tuple
 *
 * The tuple-generating function should exploit tuples returning references
 * using a form like:
 *
 * @code
 *   auto tie() const {
 *     // In order of desired lexicographical comparison
 *     return std::tie(member1, member2, ...);
 *   }
 * @endcode
 *
 */
template<typename T>
struct LexicographicComparable {
  static constexpr const T& getDerived(const LexicographicComparable& base) {
    return static_cast<const T&>(base);
  }

  constexpr bool operator == (const LexicographicComparable& other) const {
    return getDerived(*this).tie() == getDerived(other).tie();
  }

  constexpr bool operator != (const LexicographicComparable& other) const {
    return getDerived(*this).tie() != getDerived(other).tie();
  }

  constexpr bool operator < (const LexicographicComparable& other) const {
    return getDerived(*this).tie() < getDerived(other).tie();
  }

  constexpr bool operator <= (const LexicographicComparable& other) const {
    return getDerived(*this).tie() <= getDerived(other).tie();
  }

  constexpr bool operator > (const LexicographicComparable& other) const {
    return getDerived(*this).tie() > getDerived(other).tie();
  }

  constexpr bool operator >= (const LexicographicComparable& other) const {
    return getDerived(*this).tie() >= getDerived(other).tie();
  }
};

} // namespace crtp
} // namespace temple
} // namespace Scine

#endif
