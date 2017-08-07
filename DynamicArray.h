#ifndef INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_ARRAY_H
#define INCLUDE_CONSTEXPR_MAGIC_DYNAMIC_ARRAY_H

#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <utility>

#include "Containers.h"

/* TODO
 * - refactor this using Array.h as base? Loads of forwarding functions, no
 *   need for iterator and constIterator classes
 */

namespace ConstexprMagic {

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

  // NOTE: Exists solely to circumvent GCC Bug 71504
  constexpr T _defaultConstruct(size_t index __attribute__((unused))) {
    return {};
  }

  template<size_t ... Inds>
  std::initializer_list<T> _makeInitializer(
    const DynamicArray& other,
    std::index_sequence<Inds...>
  ) {
    return {
      other[Inds]...
    };
  }

public:
  /*! Delegate constructor using another DynamicArray and an index_sequence to directly
   * form the array mem-initializer with a parameter pack expansion
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

  //! Delegate std::array ctor, using same trick as copy ctor
  template<size_t ... Inds>
  constexpr DynamicArray(
    const std::array<T, nItems>& other,
    std::index_sequence<Inds...>
  ) : _items {other[Inds]...},
      _count(other.size())
  {}

  //! Construct from std::array using same trick as copy ctor
  constexpr DynamicArray(const std::array<T, nItems>& other)
    : DynamicArray(other, std::make_index_sequence<nItems>{})
  {}

  // NOTE: Exists solely to circumvent GCC Bug 71504
  template<size_t ... Inds>
  constexpr DynamicArray(std::index_sequence<Inds...>)
    : _items { _defaultConstruct(Inds) ...},
      _count {0}
  {}

  // NOTE: Exists solely to circumvent GCC Bug 71504
  constexpr DynamicArray() : DynamicArray(std::make_index_sequence<nItems>{}) {}

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
      if(_items[i] != other[i]) {
        return false;
      }
    }

    return true;
  }

  constexpr bool operator != (const DynamicArray& other) const {
    return !(*this == other);
  }

  constexpr bool operator < (const DynamicArray& other) const {
    if(_count < other._count) {
      return true;
    }

    for(unsigned i = 0; i < _count; ++i) {
      if(_items[i] < other[i]) {
        return true;
      }
    }

    return false;
  }

  constexpr bool operator > (const DynamicArray& other) const {
    return other < *this;
  }

  // Begin and end iterators
  using BaseIteratorType = std::iterator<
    std::bidirectional_iterator_tag, // iterator category
    T,                               // value_type
    unsigned,                        // difference_type
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

  using ConstBaseIteratorType = std::iterator<
    std::bidirectional_iterator_tag, // iterator category
    T,                               // value_type
    unsigned,                        // difference_type
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
      _baseRef = other._baseRef;
      _position = other._position;

      return *this;
    }

    constexpr constIterator& operator ++ () {
      _position += 1;
      return *this;
    }

    constexpr constIterator operator ++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr constIterator& operator --() {
      _position -= 1;
      return *this;
    }

    constexpr constIterator operator -- (int) {
      iterator retval = *this;
      --(*this);
      return retval;
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

  constexpr constIterator begin() const {
    return constIterator(*this, 0);
  }

  constexpr constIterator end() const {
    return constIterator(*this, _count);
  }

  constexpr operator std::array<T, nItems> () const {
    return makeArray(std::make_index_sequence<nItems>{});
  }

  constexpr std::array<T, nItems> getDynamicArray() const {
    return makeArray(std::make_index_sequence<nItems>{});
  }
};

} // namespace ConstexprMagic

#endif
