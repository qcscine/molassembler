#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_ARRAY_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_ARRAY_H

#include <cstddef>
#include <type_traits>
#include <utility>

#include "Containers.h"

/*! @file
 *
 * A constexpr fixed-maximum size managed array so that insertions and deletions
 * do not change the type signature. Principally similar to std::vector except
 * that the maximum size must be known at compile time and cannot change.
 */

namespace temple {

template<typename T, size_t nItems>
class DynamicArray {
private:
  T _items[nItems];
  size_t _count;

  template<size_t ... Inds>
  std::array<T, nItems> makeArray(std::index_sequence<Inds...>) {
    return {{
      _items[Inds]...
    }};
  }

public:
  /*!
   * Delegate constructor using another DynamicArray and an index_sequence to
   * directly form the array mem-initializer with a parameter pack expansion
   */
  template<size_t ... Inds>
  constexpr DynamicArray(const DynamicArray& other, std::index_sequence<Inds...>) 
    :_items {other[Inds]...},
     _count(other._count)
  {}

  /* Constructing from another dynamic array is tricky since we're technically
   * not allowed to edit _items in-class, so we delegate to the previous
   * constructor and directly form the mem-initializer
   */
  constexpr DynamicArray(const DynamicArray& other) 
    : DynamicArray(other, std::make_index_sequence<nItems>{})
  {}

  //! Construct from any-size array-like container using same trick as copy ctor
  template<
    template<typename, size_t> class ArrayType,
    size_t N,
    size_t ... Inds
  > constexpr DynamicArray(
    const ArrayType<T, N>& other,
    std::index_sequence<Inds...>
  ) : _items {other.at(Inds)...},
      _count(N)
  {}

  //! Construct from any size of other array-like classes
  template<
    template<typename, size_t> class ArrayType,
    size_t N
  > constexpr DynamicArray(
    const ArrayType<T, N>& other,
    std::enable_if_t<(N <= nItems)>* = 0 // Only possible for some sizes
  ) : DynamicArray(other, std::make_index_sequence<N>{}) 
  {}

  constexpr DynamicArray() : _items {}, _count(0) {}

  //! Parameter pack constructor, will work as long as the arguments are castable
  template<typename ...Args>
  constexpr DynamicArray(Args... args) 
    : _items {static_cast<T>(args)...},
      _count(sizeof...(args))
  {}

  constexpr void push_back(const T& item) {
    if(_count < nItems) {
      _items[_count] = item;
      _count += 1;
    } else {
      throw "Dynamic array is already full!";
    }
  }

  constexpr void push_back(T&& item) {
    if(_count < nItems) {
      _items[_count] = std::move(item);
      _count += 1;
    }
  }

  constexpr void pop_back() {
    if(_count > 0) {
      _count -= 1;
    }
  }

  constexpr void pop_back(const unsigned& numberToPop) {
    if(_count > numberToPop) {
      _count -= numberToPop;
    }
  }

  constexpr DynamicArray<T, nItems> splice(const unsigned& fromIndex) {
    DynamicArray<T, nItems> spliced;

    for(unsigned i = fromIndex; i < _count; ++i) {
      spliced.push_back(std::move(_items[i]));
    }

    pop_back(_count - fromIndex);

    return spliced;
  }

  constexpr bool validIndex(const unsigned& index) {
    return(index < _count);
  }

  constexpr T& operator[] (const unsigned& index) {
    return _items[index];
  }

  constexpr const T& operator[] (const unsigned& index) const {
    return _items[index];
  }

  constexpr T& at(const unsigned& index) {
    return _items[index];
  }

  constexpr const T& at(const unsigned& index) const {
    return _items[index];
  }

  constexpr T& front() {
    return _items[0];
  }

  constexpr const T& front() const {
    return _items[0];
  }

  constexpr T& back() {
    if(_count == 0) {
      return front();
    }

    return _items[_count - 1];
  }

  constexpr const T& back() const {
    if(_count == 0) {
      return front();
    }

    return _items[_count - 1];
  }

  constexpr size_t size() const {
    return _count;
  }

  constexpr void clear() {
    _count = 0;
  }

  constexpr bool operator == (const DynamicArray& other) const {
    if(_count != other._count) {
      return false;
    }

    for(unsigned i = 0; i < _count; ++i) {
      if(_items[i] != other._items[i]) {
        return false;
      }
    }

    return true;
  }

  constexpr bool operator != (const DynamicArray& other) const {
    if(_count == other._count) {
      for(unsigned i = 0; i < _count; ++i) {
        if(_items[i] != other._items[i]) {
          return true;
        }
      }
    }

    return false;
  }

  constexpr bool operator < (const DynamicArray& other) const {
    if(_count < other._count) {
      return true;
    }

    for(unsigned i = 0; i < _count; ++i) {
      if(_items[i] < other._items[i]) {
        return true;
      } else if(_items[i] > other._items[i]) {
        return false;
      }
    }

    return false;
  }

  constexpr bool operator > (const DynamicArray& other) const {
    return other < *this;
  }

  // Begin and end iterators
  using BaseIteratorType = std::iterator<
    std::random_access_iterator_tag, // iterator category
    T,                               // value_type
    int,                             // difference_type
    const T*,                        // pointer
    T&                               // reference
  >;

  class iterator : public BaseIteratorType {
  private:
    DynamicArray& _baseRef;
    unsigned _position;

  public:
    constexpr explicit iterator(
      DynamicArray& instance,
      unsigned&& initPosition
    ) : _baseRef(instance),
        _position(initPosition) 
    {}

    constexpr iterator(const iterator& other) 
      : _baseRef(other._baseRef),
        _position(other._position)
    {}

    constexpr iterator& operator = (const iterator& other) { 
      _baseRef = other._baseRef;
      _position = other._position;

      return *this;
    }

    constexpr iterator& operator ++ () {
      _position += 1;
      return *this;
    }

    constexpr iterator operator ++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr iterator& operator --() {
      _position -= 1;
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
      _position += increment;
      return *this;
    }

    constexpr iterator& operator -= (const int& increment) {
      _position -= increment;
      return *this;
    }

    constexpr int operator - (const iterator& other) const {
      return (
        static_cast<int>(_position)
        - static_cast<int>(other._position)
      );
    }

    constexpr bool operator == (const iterator& other) const {
      return (
        &_baseRef == &other._baseRef
        && _position == other._position
      );
    }

    constexpr bool operator != (const iterator& other) const {
      return !(
        *this == other
      );
    }

    constexpr typename BaseIteratorType::reference operator * () const {
      return _baseRef[_position];
    }
  };

  constexpr iterator begin() {
    return iterator(*this, 0);
  }

  constexpr iterator end() {
    return iterator(*this, _count);
  }

private:
  constexpr void _moveElementsRightUntil(
    const iterator& position
  ) {
    // Add the last element in the array onto the end
    push_back(
      back()
    );

    // Move everything up to and including the lower bound back
    auto rightIter = end();
    rightIter -= 2;
    auto leftIter = rightIter - 1;

    while(rightIter != position) {
      *rightIter = std::move(*leftIter);

      --leftIter;
      --rightIter;
    }
  }

public:
  constexpr void insertAt(
    iterator insertPositionIter,
    const T& item
  ) {
    // In case there is no object larger than the one being inserted, add to end
    if(insertPositionIter == end()) {
      push_back(item);
    } else {
      _moveElementsRightUntil(insertPositionIter);

      // Copy in the item
      *insertPositionIter = item;
    }
  }

  //! Inserts an element at a specified position O(N)
  constexpr void insertAt(
    iterator insertPositionIter,
    T&& item
  ) {
    // In case there is no object larger than the one being inserted, add to end
    if(insertPositionIter == end()) {
      push_back(item);
    } else {
      _moveElementsRightUntil(insertPositionIter);

      // Move in the item
      *insertPositionIter = std::move(item);
    }
  }

  constexpr void removeAt(iterator insertPositionIter) {
    if(insertPositionIter == end()) {
      throw "Cannot remove item at end iterator!";
    }

    // Rename for clarity
    auto& leftIter = insertPositionIter;
    auto rightIter = insertPositionIter + 1;

    while(rightIter != end()) {
      *leftIter = std::move(*rightIter);

      ++leftIter;
      ++rightIter;
    }

    pop_back();
  }

  using ConstBaseIteratorType = std::iterator<
    std::random_access_iterator_tag, // iterator category
    T,                               // value_type
    int,                             // difference_type
    const T*,                        // pointer
    const T&                         // reference
  >;
  
  class constIterator : public ConstBaseIteratorType {
  private:
    const DynamicArray& _baseRef;
    unsigned _position;

  public:
    constexpr explicit constIterator(
      const DynamicArray& instance,
      unsigned&& initPosition
    ) : _baseRef(instance),
        _position(initPosition) 
    {}

    constexpr constIterator(const constIterator& other) 
      : _baseRef(other._baseRef),
        _position(other._position)
    {}

    constexpr constIterator& operator = (const constIterator& other) { 
      if(_baseRef != other._baseRef) {
        throw "Trying to assign constIterator to other base DynamicArray!";
      }

      _position = other._position;

      return *this;
    }

    constexpr constIterator& operator ++ () {
      _position += 1;
      return *this;
    }

    constexpr constIterator operator ++ (int) {
      constIterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr constIterator& operator --() {
      _position -= 1;
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
      _position += increment;
      return *this;
    }

    constexpr constIterator& operator -= (const int& increment) {
      _position -= increment;
      return *this;
    }

    constexpr int operator - (const constIterator& other) const {
      return (
        static_cast<int>(_position)
        - static_cast<int>(other._position)
      );
    }

    constexpr bool operator == (const constIterator& other) const {
      return (
        &_baseRef == &other._baseRef
        && _position == other._position
      );
    }

    constexpr bool operator != (const constIterator& other) const {
      return !(
        *this == other
      );
    }

    constexpr typename ConstBaseIteratorType::reference operator * () const {
      return _baseRef[_position];
    }
  };

  //! Type alias for compatibility with STL algorithms
  using const_iterator = constIterator;

  constexpr constIterator begin() const {
    return constIterator(*this, 0);
  }

  constexpr constIterator end() const {
    return constIterator(*this, _count);
  }

  constexpr operator std::array<T, nItems> () const {
    return makeArray(std::make_index_sequence<nItems>{});
  }

  constexpr void copyIn(const DynamicArray<T, nItems>& other) {
    if(other.size() + size() > nItems) {
      throw "DynamicArray to be copied in has too many elements to fit!";
    }

    for(auto it = other.begin(); it != other.end(); ++it) {
      push_back(*it);
    }
  }
};

/*!
 * Groups data by equality.
 */
template<
  template<typename, size_t> class ArrayType,
  typename T,
  size_t N,
  class BinaryFunction
> constexpr DynamicArray<
  DynamicArray<T, N>,
  N
> groupByEquality(
  const ArrayType<T, N>& data,
  BinaryFunction&& equalityComparator
) {
  // Maximal dimensions if all equal is 1xN, if all unequal Nx1
  DynamicArray<
    DynamicArray<T, N>,
    N
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
        DynamicArray<T, N> {*iter}
      );
    }
  }

  return groups;
}

template<typename T, size_t N>
DynamicArray<T, N> merge(
  const DynamicArray<T, N>& a,
  const DynamicArray<T, N>& b
) {
  if(a.size() + b.size() > N) {
    throw "DynamicArrays to be merged have too many elements!";
  }

  DynamicArray<T, N> merged {a};
  merged.copyIn(b);
  return merged;
}

} // namespace temple

#endif
