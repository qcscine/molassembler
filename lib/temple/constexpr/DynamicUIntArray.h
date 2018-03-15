#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_UINT_ARRAY_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_UINT_ARRAY_H

#include "Math.h"
#include "Array.h"
#include "DynamicArray.h"

/*! @file
 *
 * Implements a dynamic vector-like class limited to small amounts of the
 * numbers 0-9 that has particularly fast comparison operators.
 */

namespace temple {

/*!
 * A dynamic vector-like class for small amounts of the numbers 0-9.
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
 * This has the effect, however, that access and moficiation are slower, and
 * some methods must return proxy object instances for array-like syntax.
 */
template<typename UnsignedType>
class DynamicUIntArray {
public:
  // Size of the array
  static constexpr unsigned N = Math::floor(
    Math::log10(
      static_cast<double>(
        std::numeric_limits<UnsignedType>::max()
      )
    )
  );

private:
  UnsignedType _data;
  unsigned _count;

  template<size_t ... Inds>
  DynamicArray<unsigned, N> _toDynamicArray(std::index_sequence<Inds...>) const {
    return {{
      at(Inds)...
    }};
  }

public:

  /* Constructors */
  constexpr DynamicUIntArray() : _data(0), _count(0) {}

  template<size_t size>
  constexpr explicit DynamicUIntArray(const Array<unsigned, size>& values) 
    : _data(0), 
      _count(static_cast<unsigned>(size)) 
  {
    static_assert(
      size <= N,
      "Size of initializing array must fit within DynamicUIntArray!"
    );

    for(unsigned i = 0; i < size; ++i) {
      _data += Math::pow(10, i) * values.at(i);
    }
  }

  template<typename ... Args>
  constexpr DynamicUIntArray(Args ... args) 
  : _data(0), 
    _count(sizeof...(args)) 
  {
    static_assert(
      sizeof...(args) <= N,
      "Variadic parameter pack must fit within DynamicUIntArray!"
    );

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
    DynamicUIntArray* _basePtr;

  public:
    constexpr ModificationProxy(
      unsigned digit,
      DynamicUIntArray* basePtr
    ) : _digit(digit), _basePtr(basePtr) 
    {
      if(digit >= basePtr->_count) {
        throw "UInt non-const Array access out of bounds";
      }
    }

    constexpr void operator = (unsigned assignValue) {
      if(assignValue > 9) {
        throw "Trying to assign a value larger than 9 in DynamicUIntArray";
      }

      _basePtr->_data += Math::pow(10, _digit) * (
        static_cast<int>(assignValue) // New value
        - const_cast<const DynamicUIntArray*>(_basePtr)->at(_digit) // Existing value
      );
    }

    constexpr operator unsigned() const {
      return const_cast<const DynamicUIntArray*>(_basePtr)->at(_digit);
    }
  };

  constexpr ModificationProxy at(const unsigned& i) {
    return {
      i,
      this
    };
  }

  constexpr void push_back(const unsigned& i) {
    _data.at(_count) = i;
    ++_count;
  }

  constexpr void pop_back() {
    _data.at(_count - 1) = 0;
    --_count;
  }

  /* Information */

  constexpr unsigned at(const unsigned& i) const {
    if(i >= _count) {
      throw "DynamicUIntArray access out of bounds";
    }

    return (_data / Math::pow(10, i)) % 10;
  }

  constexpr unsigned front() const {
    return at(0);
  }

  constexpr unsigned back() const {
    return at(_count - 1);
  }

  constexpr unsigned size() const {
    return _count;
  }


  /* Operators */
  constexpr bool operator == (const DynamicUIntArray& other) const {
    return _data == other._data;
  }

  constexpr bool operator != (const DynamicUIntArray& other) const {
    return _data != other._data;
  }

  constexpr bool operator < (const DynamicUIntArray& other) const {
    return _data > other._data;
  }

  constexpr bool operator > (const DynamicUIntArray& other) const {
    return _data < other._data;
  }


  /* Iterators */

  using BaseIteratorType = std::iterator<
    std::random_access_iterator_tag, // iterator category
    unsigned,                        // value_type
    int,                             // difference_type
    unsigned,                        // pointer
    ModificationProxy                // reference
  >;

  class iterator : public BaseIteratorType {
  private:
    unsigned _digit;
    DynamicUIntArray* _basePtr;

  public:
    constexpr iterator(
      unsigned digit,
      DynamicUIntArray* basePtr
    ) : _digit(digit), _basePtr(basePtr) 
    {
      if(digit > basePtr->_count) {
        throw "Initialization of iterator out of bounds of DynamicUIntArray parent";
      }
    }

    constexpr iterator(const iterator& other)
      : _digit(other._digit),
        _basePtr(other._basePtr)
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
        _basePtr == other._basePtr
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
    return iterator {0, this};
  }

  constexpr iterator end() {
    return iterator {_count, this};
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
    const DynamicUIntArray* const _basePtr;

  public:
    constexpr constIterator(
      unsigned digit,
      const DynamicUIntArray* const basePtr
    ) : _digit(digit),
        _basePtr(basePtr) 
    {
      if(digit > basePtr->_count) {
        throw "Initialization of const iterator out of bounds of DynamicUIntArray parent";
      }
    }

    constexpr constIterator(const constIterator& other) 
      : _digit(other._digit),
        _basePtr(other._basePtr)
    {}

    constexpr constIterator& operator = (const constIterator& other) { 
      if(_basePtr != other._basePtr) {
        throw "Trying to assign constIterator to other base DynamicArray!";
      }

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
        _basePtr == other._basePtr
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
    return constIterator {0, this};
  }

  constexpr constIterator end() const {
    return constIterator {_count, this};
  }

  DynamicArray<unsigned, N> toDynamicArray() const {
    return _toDynamicArray(std::make_index_sequence<N>());
  }
};

/*!
 * Groups data by equality.
 */
template<
  typename UnsignedType,
  class BinaryFunction
> constexpr auto groupByEquality(
  const DynamicUIntArray<UnsignedType>& data,
  BinaryFunction&& equalityComparator
) {
  constexpr size_t maxSize = DynamicUIntArray<UnsignedType>::N;

  // Maximal dimensions if all equal is 1xN, if all unequal Nx1
  DynamicArray<
    DynamicArray<unsigned, maxSize>,
    maxSize
  > groups;

  for(auto iter = data.begin(); iter != data.end(); ++iter) {
    bool foundEqual = false;

    for(auto& group : groups) {
      if(equalityComparator(*iter, *group.begin())) {
        group.push_back(*iter);
        foundEqual = true;
        break;
      }
    }

    if(!foundEqual) {
      groups.push_back(
        DynamicArray<unsigned, maxSize> {*iter}
      );
    }
  }

  return groups;
}

} // namespace temple

#endif
