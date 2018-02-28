#ifndef INCLUDE_CONSTEXPR_MAGIC_OPTIONAL_H
#define INCLUDE_CONSTEXPR_MAGIC_OPTIONAL_H

#include <utility>

/*! @file
 *
 * Includes an optional type implementation.
 */

namespace constable {

/*!
 * A constexpr option type much like std::optional with the limitation that T 
 * must be default-constructible and equality and less-than comparison operators
 * must be defined for the underlying type.
 */
template<typename T>
class Optional {
private:
  T _value;
  bool _hasValue;

public:
  constexpr Optional() : _value(T {}), _hasValue(false) {}
  constexpr Optional(T value) : _value(std::move(value)), _hasValue(true) {}

  constexpr operator bool () const {
    return _hasValue;
  }

  constexpr bool hasValue() const {
    return _hasValue;
  }

  constexpr T value() const {
    return _value;
  }

  constexpr T valueOr(const T& alternative) const {
    if(_hasValue) {
      return _value;
    }

    return alternative;
  }

  constexpr void operator = (T assignment) {
    _value = assignment;
  }

  constexpr bool operator == (const Optional& other) const {
    if(!_hasValue && !other._hasValue) {
      return true;
    }

    if(_hasValue && other._hasValue) {
      return _value == other._value;
    }

    return false;
  }

  constexpr bool operator != (const Optional& other) const {
    return !(*this == other);
  }

  constexpr bool operator < (const Optional& other) const {
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

  constexpr bool operator > (const Optional& other) const {
    return (other < *this);
  }
};

} // namespace constable

#endif
