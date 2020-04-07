/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Basic optional type
 */

#ifndef INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_OPTIONAL_H
#define INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_OPTIONAL_H

#include "molassembler/Temple/Preprocessor.h"

#include <utility>

namespace Scine {
namespace temple {

/*! @brief An Option monadic type
 *
 * A constexpr option type much like std::optional with the limitation that T
 * must be:
 * - DefaultConstructible
 * - LessThanComparable
 * - EqualityComparable
 *
 * (see the C++ standard definitions)
 *
 * @warning T may not be a reference or const-qualified type. Those are
 * impossible to implement constexpr in C++14 as far as I can tell.
 */
template<typename T>
class Optional {
public:
//!@name Constructors
//!@{
  /*! @brief Default constructor
   *
   * None value. Default constructs a contained type.
   */
  constexpr Optional() {
    static_assert(
      std::is_same<T, std::decay_t<T>>::value,
      "T must not be a reference or const-qualified type"
    );
  }

  //! Value constructor
  constexpr explicit Optional(T value) : value_(std::move(value)), hasValue_(true) {
    static_assert(
      std::is_same<T, std::decay_t<T>>::value,
      "T must not be a reference or const-qualified type"
    );
  }
//!@}

//!@name Information
//!@{
  /*! @brief Returns whether the optional contains a value
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr bool hasValue() const {
    return hasValue_;
  }

  /*! @brief Returns the contained value unchecked
   *
   * @complexity{@math{\Theta(1)}}
   * @warning If @p hasValue is false, this is UB.
   */
  PURITY_WEAK constexpr T value() const {
    return value_;
  }

  /*! @brief Monadic bind with function of signature T -> U
   *
   * @tparam UnaryFunction: Function of signature T -> U
   *
   * @returns Optional<U>
   */
  template<class UnaryFunction>
  constexpr auto map(UnaryFunction&& function) const {
    // Function has signature T -> U
    using U = decltype(function(value_));

    if(hasValue_) {
      return Optional<U> {
        function(value_)
      };
    }

    return Optional<U> {};
  }

  /*! @brief Monadic bind with function of signature T -> Optional<U>
   *
   * @tparam UnaryFunction: Function of signature T -> Optional<U>
   *
   * @returns Optional<U>
   */
  template<class UnaryFunction>
  constexpr auto flatMap(UnaryFunction&& function) const {
    // Function has signature T -> Optional<U>
    using OptionalU = decltype(function(value_));

    if(hasValue_) {
      return OptionalU {function(value_)};
    }

    return OptionalU {};
  }

  /*! @brief Returns a value if initialized, and another if not
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr T valueOr(const T& alternative) const {
    if(hasValue_) {
      return value_;
    }

    return alternative;
  }
//!@}

//!@name Operators
//!@{
  //! Assignment from T
  constexpr Optional& operator = (T assignment) {
    value_ = assignment;
    hasValue_ = true;

    return *this;
  }

  //! Convert-to-bool operator
  PURITY_WEAK constexpr operator bool () const {
    return hasValue_;
  }

  //! Compares on basis of contained value. Nones do compare equal
  PURITY_WEAK constexpr bool operator == (const Optional& other) const {
    if(!hasValue_ && !other.hasValue_) {
      return true;
    }

    if(hasValue_ && other.hasValue_) {
      return value_ == other.value_;
    }

    return false;
  }

  //! Compares on basis of contained value. Nones do compare equal
  PURITY_WEAK constexpr bool operator != (const Optional& other) const {
    return !(*this == other);
  }

  //! Lexicographical-like comparison
  PURITY_WEAK constexpr bool operator < (const Optional& other) const {
    // If neither has a value, they are equal
    if(!hasValue_ && !other.hasValue_) {
      return false;
    }

    // If we do not have a value, but the other does, we are smaller
    if(!hasValue_ && other.hasValue_) {
      return true;
    }

    // If we have a value, but the other doesn't, we are bigger
    if(hasValue_ && !other.hasValue_) {
      return false;
    }

    // Remaining case: both have values
    return (
      value_ < other.value_
    );
  }

  //! Lexicographical-like comparison
  PURITY_WEAK constexpr bool operator > (const Optional& other) const {
    return (other < *this);
  }
//!@}

private:
//!@name State
//!@{
  T value_ = T {};
  bool hasValue_ = false;
//!@}
};

template<typename T>
class Optional<T&> {
public:
//!@name Constructors
//!@{
  /*! @brief Default constructor
   *
   * None value. Default constructs a contained type.
   */
  constexpr Optional() : ref_(dummy_) {}

  //! Value constructor
  constexpr explicit Optional(T& value) : ref_(value), hasValue_(true) {}
//!@}

//!@name Information
//!@{
  /*! @brief Returns whether the optional contains a value
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr bool hasValue() const {
    return hasValue_;
  }

  /*! @brief Returns the contained value unchecked
   *
   * @complexity{@math{\Theta(1)}}
   * @warning If @p hasValue is false, this is UB.
   */
  PURITY_WEAK constexpr T& value() const {
    return ref_;
  }

  /*! @brief Monadic bind with function of signature T -> U
   *
   * @tparam UnaryFunction: Function of signature T -> U
   *
   * @returns Optional<U>
   */
  template<class UnaryFunction>
  constexpr auto map(UnaryFunction&& function) const {
    // Function has signature T -> U
    using U = decltype(function(ref_));

    if(hasValue_) {
      return Optional<U> {
        function(ref_)
      };
    }

    return Optional<U> {};
  }

  /*! @brief Monadic bind with function of signature T -> Optional<U>
   *
   * @tparam UnaryFunction: Function of signature T -> Optional<U>
   *
   * @returns Optional<U>
   */
  template<class UnaryFunction>
  constexpr auto flatMap(UnaryFunction&& function) const {
    // Function has signature T -> Optional<U>
    using OptionalU = decltype(function(ref_));

    if(hasValue_) {
      return OptionalU {function(ref_)};
    }

    return OptionalU {};
  }

  /*! @brief Returns a value if initialized, and another if not
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr T valueOr(const T& alternative) const {
    if(hasValue_) {
      return ref_;
    }

    return alternative;
  }
//!@}

//!@name Operators
//!@{
  //! Convert-to-bool operator
  PURITY_WEAK constexpr operator bool () const {
    return hasValue_;
  }

  //! Compares on basis of contained value. Nones do compare equal
  PURITY_WEAK constexpr bool operator == (const Optional& other) const {
    if(!hasValue_ && !other.hasValue_) {
      return true;
    }

    if(hasValue_ && other.hasValue_) {
      return ref_ == other.ref_;
    }

    return false;
  }

  //! Compares on basis of contained value. Nones do compare equal
  PURITY_WEAK constexpr bool operator != (const Optional& other) const {
    return !(*this == other);
  }

  //! Lexicographical-like comparison
  PURITY_WEAK constexpr bool operator < (const Optional& other) const {
    // If neither has a value, they are equal
    if(!hasValue_ && !other.hasValue_) {
      return false;
    }

    // If we do not have a value, but the other does, we are smaller
    if(!hasValue_ && other.hasValue_) {
      return true;
    }

    // If we have a value, but the other doesn't, we are bigger
    if(hasValue_ && !other.hasValue_) {
      return false;
    }

    // Remaining case: both have values
    return (
      ref_ < other.ref_
    );
  }

  //! Lexicographical-like comparison
  PURITY_WEAK constexpr bool operator > (const Optional& other) const {
    return (other < *this);
  }
//!@}

private:
//!@name State
//!@{
  T dummy_ = T{};
  T& ref_;
  bool hasValue_ = false;
//!@}
};

} // namespace temple
} // namespace Scine

#endif
