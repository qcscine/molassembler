/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Basic optional type
 */

#ifndef INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_OPTIONAL_H
#define INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_OPTIONAL_H

#include "temple/Preprocessor.h"

#include <utility>

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
  constexpr explicit Optional(T value) : _value(std::move(value)), _hasValue(true) {
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
    return _hasValue;
  }

  /*! @brief Returns the contained value unchecked
   *
   * @complexity{@math{\Theta(1)}}
   * @warning If @p hasValue is false, this is UB.
   */
  PURITY_WEAK constexpr T value() const {
    return _value;
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
    using U = decltype(function(_value));

    if(_hasValue) {
      return Optional<U> {
        function(_value)
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
    using OptionalU = decltype(function(_value));

    if(_hasValue) {
      return OptionalU {function(_value)};
    }

    return OptionalU {};
  }

  /*! @brief Returns a value if initialized, and another if not
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr T valueOr(const T& alternative) const {
    if(_hasValue) {
      return _value;
    }

    return alternative;
  }
//!@}

//!@name Operators
//!@{
  //! Assignment from T
  constexpr Optional& operator = (T assignment) {
    _value = assignment;
    _hasValue = true;

    return *this;
  }

  //! Convert-to-bool operator
  PURITY_WEAK constexpr operator bool () const {
    return _hasValue;
  }

  //! Compares on basis of contained value. Nones do compare equal
  PURITY_WEAK constexpr bool operator == (const Optional& other) const {
    if(!_hasValue && !other._hasValue) {
      return true;
    }

    if(_hasValue && other._hasValue) {
      return _value == other._value;
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
    if(!_hasValue && !other._hasValue) {
      return false;
    }

    // If we do not have a value, but the other does, we are smaller
    if(!_hasValue && other._hasValue) {
      return true;
    }

    // If we have a value, but the other doesn't, we are bigger
    if(_hasValue && !other._hasValue) {
      return false;
    }

    // Remaining case: both have values
    return (
      _value < other._value
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
  T _value = T {};
  bool _hasValue = false;
//!@}
};

} // namespace temple

#endif
