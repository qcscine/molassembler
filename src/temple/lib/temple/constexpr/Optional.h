#ifndef INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_OPTIONAL_H
#define INLCUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_OPTIONAL_H

#include "temple/Preprocessor.h"

#include <utility>

/*! @file
 *
 * Includes an optional type implementation.
 */

namespace temple {

/*!
 * A constexpr option type much like std::optional with the limitation that T
 * must be default-constructible and equality and less-than comparison operators
 * must be defined for the underlying type.
 *
 * @warning T must not be reference or const-qualified type
 */
template<typename T>
class Optional {
public:
//!@name Constructors
//!@{
  //! Default constructor
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
  //! Returns whether the optional contains a value
  constexpr bool hasValue() const PURITY_WEAK {
    return _hasValue;
  }

  //! Returns a value unchecked
  constexpr T value() const PURITY_WEAK {
    return _value;
  }

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

  template<class UnaryFunction>
  constexpr auto flatMap(UnaryFunction&& function) const {
    // Function has signature T -> Optional<U>
    using OptionalU = decltype(function(_value));

    if(_hasValue) {
      return OptionalU {function(_value)};
    }

    return OptionalU {};
  }

  //! Returns a value if initialized, and another if not
  constexpr T valueOr(const T& alternative) const PURITY_WEAK {
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
  constexpr operator bool () const PURITY_WEAK {
    return _hasValue;
  }

  //! Compares on basis of contained value. Nones do compare equal
  constexpr bool operator == (const Optional& other) const PURITY_WEAK {
    if(!_hasValue && !other._hasValue) {
      return true;
    }

    if(_hasValue && other._hasValue) {
      return _value == other._value;
    }

    return false;
  }

  //! Compares on basis of contained value. Nones do compare equal
  constexpr bool operator != (const Optional& other) const PURITY_WEAK {
    return !(*this == other);
  }

  //! Lexicographical-like comparison
  constexpr bool operator < (const Optional& other) const PURITY_WEAK {
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
  constexpr bool operator > (const Optional& other) const PURITY_WEAK {
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
