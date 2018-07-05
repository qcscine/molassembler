#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_UINT_ARRAY_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_UINT_ARRAY_H

#include "Math.h"
#include "Array.h"

/*! @file
 *
 * Implements a fixed-size array-like class limited to small amounts of the
 * numbers 0-9 that has particularly fast comparison operators.
 */

namespace temple {

/*!
 * A fixed array-like class for small amounts of the numbers 0-9.
 *
 * The following amount of numbers can be stored for the template unsigned
 * integer types:
 *
 * - 8 bit -> 2
 * - 16 bit -> 5
 * - 32 bit -> 9
 * - 64 bit -> 19
 *
 * The numbers 0-9 are stored as digits in an unsigned decimal number. So
 * although the class represents several unsigned integers, all operators are
 * evaluated on a single number, which is O(1).
 *
 * The tradeoff is that access and modification are slower.
 */
template<typename UnsignedType>
class UIntArray {
private:
  UnsignedType _data;

public:
  // Size of the array
  static constexpr unsigned N = Math::floor(
    Math::log10(
      static_cast<double>(
        std::numeric_limits<UnsignedType>::max()
      )
    )
  );


  /* Constructors */
  constexpr UIntArray() : _data(0) {}

  template<size_t size>
  constexpr UIntArray(const Array<unsigned, size>& values) : _data(0) {
    static_assert(
      size <= N,
      "Size of initializing array must fit within UIntArray!"
    );

    for(unsigned i = 0; i < size; ++i) {
      _data += Math::pow(10, i) * values.at(i);
    }
  }

  template<typename ... Args>
  constexpr UIntArray(Args ... args) : _data(0) {
    Array<UnsignedType, sizeof...(args)> values {{
      static_cast<UnsignedType>(args)...
    }};

    for(unsigned i = 0; i < sizeof...(args); ++i) {
      _data += Math::pow(10, i) * values.at(i);
    }
  }


  /* Modification */
  class ModificationProxy {
  private:
    unsigned _digit;
    UIntArray* _basePtr;

  public:
    constexpr ModificationProxy(
      unsigned digit,
      UIntArray* basePtr
    ) : _digit(digit), _basePtr(basePtr)
    {
      if(digit >= N) {
        throw "UInt non-const Array access out of bounds";
      }
    }

    constexpr void operator = (unsigned assignValue) {
      if(assignValue > 9) {
        throw "Trying to assign a value larger than 9 in UIntArray";
      }

      _basePtr->_data += Math::pow(10, _digit) * (
        static_cast<int>(assignValue) // New value
        - const_cast<const UIntArray*>(_basePtr)->at(_digit) // Existing value
      );
    }

    constexpr operator unsigned() const {
      return const_cast<const UIntArray*>(_basePtr)->at(_digit);
    }
  };

  constexpr ModificationProxy at(const unsigned& i) {
    return {
      i,
      this
    };
  }

  /* Information */

  constexpr unsigned at(const unsigned& i) const {
    if(i >= N) {
      throw "UIntArray access out of bounds";
    }

    return (_data / Math::pow(10, i)) % 10;
  }

  constexpr unsigned front() const {
    return at(0);
  }

  constexpr unsigned back() const {
    return at(N - 1);
  }


  /* Operators */
  constexpr bool operator == (const UIntArray& other) const {
    return _data == other._data;
  }

  constexpr bool operator != (const UIntArray& other) const {
    return _data != other._data;
  }

  constexpr bool operator < (const UIntArray& other) const {
    return _data > other._data;
  }

  constexpr bool operator > (const UIntArray& other) const {
    return _data < other._data;
  }


  /* Iterators */

  using BaseIteratorType = std::iterator<
    std::random_access_iterator_tag, // iterator category
    unsigned,                        // value_type
    int,                             // difference_type
    unsigned,                        // pointer
    ModificationProxy&               // reference
  >;

  class iterator : public BaseIteratorType {
  private:
    unsigned _digit;
    UIntArray* _basePtr;

  public:
    constexpr iterator(
      unsigned digit,
      UIntArray* basePtr
    ) : _digit(digit), _basePtr(basePtr)
    {
      if(digit > N) {
        throw "Initialization of iterator out of bounds of UIntArray parent";
      }
    }

    constexpr iterator(const iterator& other)
      : _basePtr(other._basePtr),
        _digit(other._digit)
    {}

    constexpr iterator& operator = (const iterator& other) {
      _basePtr = other._basePtr;
      _digit = other._digit;
    }

    constexpr iterator& operator ++ () {
      _digit += 1;
      return *this;
    }

    constexpr iterator operator ++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr iterator& operator --() {
      _digit -= 1;
      return *this;
    }

    constexpr iterator operator -- (int) {
      iterator retval = *this;
      --(*this);
      return retval;
    }

    constexpr iterator operator + (const int& increment) {
      iterator retval = *this;
      retval += increment;
      return retval;
    }

    constexpr iterator operator - (const int& increment) {
      iterator retval = *this;
      retval -= increment;
      return retval;
    }

    constexpr iterator& operator += (const int& increment) {
      _digit += increment;
      return *this;
    }

    constexpr iterator& operator -= (const int& increment) {
      _digit -= increment;
      return *this;
    }

    constexpr int operator - (const iterator& other) const {
      return (
        static_cast<int>(_digit)
        - static_cast<int>(other._digit)
      );
    }

    constexpr bool operator == (const iterator& other) const {
      return (
        &_basePtr == &other._basePtr
        && _digit == other._digit
      );
    }

    constexpr bool operator != (const iterator& other) const {
      return !(
        *this == other
      );
    }

    constexpr typename BaseIteratorType::reference operator * () const {
      return _basePtr->at(_digit);
    }
  };

  constexpr iterator begin() {
    return iterator {0, *this};
  }

  constexpr iterator end() {
    return iterator {N, *this};
  }

  using ConstBaseIteratorType = std::iterator<
    std::random_access_iterator_tag, // iterator category
    unsigned,                        // value_type
    int,                             // difference_type
    unsigned,                        // pointer
    unsigned                         // reference
  >;

  class constIterator : public ConstBaseIteratorType {
  private:
    unsigned _digit;
    UIntArray* const _basePtr;

  public:
    constexpr constIterator(
      unsigned digit,
      UIntArray* const basePtr
    ) : _digit(digit), _basePtr(basePtr)
    {
      if(digit > N) {
        throw "Initialization of const iterator out of bounds of UIntArray parent";
      }
    }

    constexpr constIterator(const constIterator& other)
      : _basePtr(other._basePtr),
        _digit(other._digit)
    {}

    constexpr constIterator& operator = (const constIterator& other) {
      _basePtr = other._basePtr;
      _digit = other._digit;

      return *this;
    }

    constexpr constIterator& operator ++ () {
      _digit += 1;
      return *this;
    }

    constexpr constIterator operator ++ (int) {
      constIterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr constIterator& operator --() {
      _digit -= 1;
      return *this;
    }

    constexpr constIterator operator -- (int) {
      constIterator retval = *this;
      --(*this);
      return retval;
    }

    constexpr constIterator operator + (const int& increment) {
      constIterator retval = *this;
      retval += increment;
      return retval;
    }

    constexpr constIterator operator - (const int& increment) {
      constIterator retval = *this;
      retval -= increment;
      return retval;
    }

    constexpr constIterator& operator += (const int& increment) {
      _digit += increment;
      return *this;
    }

    constexpr constIterator& operator -= (const int& increment) {
      _digit -= increment;
      return *this;
    }

    constexpr int operator - (const constIterator& other) const {
      return (
        static_cast<int>(_digit)
        - static_cast<int>(other._digit)
      );
    }

    constexpr bool operator == (const constIterator& other) const {
      return (
        &_basePtr == &other._basePtr
        && _digit == other._digit
      );
    }

    constexpr bool operator != (const constIterator& other) const {
      return !(
        *this == other
      );
    }

    constexpr typename ConstBaseIteratorType::reference operator * () const {
      return _basePtr->at(_digit);
    }
  };

  //! Type alias for compatibility with STL algorithms
  using const_iterator = constIterator;

  constexpr constIterator begin() const {
    return constIterator(0, this);
  }

  constexpr constIterator end() const {
    return constIterator(N, this);
  }
};

} // namespace temple

#endif
