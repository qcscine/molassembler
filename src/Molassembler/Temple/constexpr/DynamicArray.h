/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief std::vector-like class (but max size is size allocated)
 *
 * A constexpr fixed-maximum size managed array so that insertions and deletions
 * do not change the type signature. Principally similar to std::vector except
 * that the maximum size must be known at compile time and cannot change.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_DYNAMIC_ARRAY_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_DYNAMIC_ARRAY_H

#include "Molassembler/Temple/Preprocessor.h"
#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace Scine {
namespace Molassembler {
namespace Temple {

template<typename T, std::size_t nItems>
class DynamicArray {
public:
  // Forward-declare the iterator types
  class iterator;
  class const_iterator;

//!@name Constructors
//!@{
  //! Default constructor
  constexpr DynamicArray() : items_ {} {}

  /*! @brief Helper to the copy constructor
   *
   * Delegate copy constructor directly forms the array mem-initializer with a
   * parameter pack expansion
   */
  template<std::size_t ... Inds>
  constexpr DynamicArray(
    const DynamicArray& other,
    std::index_sequence<Inds...> /* inds */
  ) :items_ {other[Inds]...},
     count_(other.count_)
  {}

  /*! @brief Copy constructor
   *
   * Constructing from another dynamic array is tricky since we're technically
   * not allowed to edit items_ in-class, so we delegate to the previous
   * constructor and directly form the mem-initializer
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr DynamicArray(const DynamicArray& other)
    : DynamicArray(other, std::make_index_sequence<nItems>{}) {}

  /* @brief Helper to the move constructor
   *
   * Delegate move constructor to directly form the array mem-initializer with
   * a parameter pack expansion
   */
  template<std::size_t ... Inds>
  constexpr DynamicArray(
    DynamicArray&& other,
    std::index_sequence<Inds...> /* inds */
  ) :items_ {std::move(other[Inds])...},
     count_(other.count_)
  {}

  /*! @brief Move constructor
   *
   * Constructed using the same technique as the copy constructor using a
   * delegate.
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr DynamicArray(DynamicArray&& other) noexcept
    : DynamicArray(std::move(other), std::make_index_sequence<nItems>{}) {}

  /*! @brief Copy assignment operator
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr DynamicArray& operator = (const DynamicArray& other) {
    for(unsigned i = 0; i < other.count_; ++i) {
      items_[i] = other.items_[i];
    }
    count_ = other.count_;
    return *this;
  }

  /*! @brief Move assignment operator
   *
   * Moves unfortunately aren't possible for the underlying array as in a
   * pointer swap, but instead must be performed with N moves of the stored
   * types.
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr DynamicArray& operator = (DynamicArray&& other) noexcept {
    for(unsigned i = 0; i < other.count_; ++i) {
      items_[i] = std::move(other.items_[i]);
    }
    count_ = other.count_;
    return *this;
  }

  ~DynamicArray() = default;
//!@}


//!@name Converting constructors
//!@{
  //! Construct from any-size array-like container using same trick as copy ctor
  template<
    template<typename, std::size_t> class ArrayType,
    std::size_t N,
    std::size_t ... Inds
  > constexpr DynamicArray(
    const ArrayType<T, N>& other,
    std::index_sequence<Inds...> /* inds */
  ) : items_ {other.at(Inds)...},
      count_(N)
  {}

  //! Construct from any size of other array-like classes
  template<
    template<typename, std::size_t> class ArrayType,
    std::size_t N
  > constexpr DynamicArray(
    const ArrayType<T, N>& other,
    std::enable_if_t<(N <= nItems)>* /* f */ = 0 // Only possible for some sizes
  ) : DynamicArray(other, std::make_index_sequence<N>{})
  {}

  //! Parameter pack constructor, will work as long as the arguments are castable
  template<typename ...Args>
  constexpr DynamicArray(Args... args)
    : items_ {static_cast<T>(args)...},
      count_(sizeof...(args))
  {}
//!@}

//!@name Modification
//!@{
  /*! @brief Adds an element to the dynamic array
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @note This does not behave like std::vector's push_back since the maximum
   * number of elements cannot be exceeded and there is never any reallocation.
   */
  constexpr void push_back(const T& item) {
    if(count_ < nItems) {
      items_[count_] = item;
      count_ += 1;
    } else {
      throw "Dynamic array is already full!";
    }
  }

  //! @overload
  constexpr void push_back(T&& item) {
    if(count_ < nItems) {
      items_[count_] = std::move(item);
      count_ += 1;
    }
  }

  /*! @brief Removes the last element
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr void pop_back() {
    if(count_ > 0) {
      count_ -= 1;
    }
  }

  /*! @brief Removes N elements from the back
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr void pop_back(const unsigned numberToPop) {
    if(count_ > numberToPop) {
      count_ -= numberToPop;
    }
  }

  /*! @brief Moves elements to a new dynamic array
   *
   * Moves items starting at a particular index to a new dynamic array.
   * This dynamic array then has fewer elements.
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr DynamicArray<T, nItems> splice(const unsigned fromIndex) {
    DynamicArray<T, nItems> spliced;

    for(unsigned i = fromIndex; i < count_; ++i) {
      spliced.push_back(std::move(items_[i]));
    }

    pop_back(count_ - fromIndex);

    return spliced;
  }
//!@}

//!@name Information
//!@{
  /*! @brief Checks whether an index is a valid accessor to underlying data
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr bool validIndex(const unsigned index) const noexcept {
    return (index < count_);
  }

  /*! @brief Returns the number of contained elements
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr std::size_t size() const noexcept {
    return count_;
  }
//!@}

//!@name Element access
//!@{
  /*! @brief Safe access to underlying data
   *
   * Possibility for UB is not "allowed" in constexpr functions, so bounds
   * checking is performed.
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr T& operator[] (const unsigned index) noexcept {
    // Defined behavior instead of UB
    if(!validIndex(index)) {
      return back();
    }

    return items_[index];
  }

  //! @overload
  PURITY_WEAK constexpr const T& operator[] (const unsigned index) const noexcept {
    if(!validIndex(index)) {
      return back();
    }

    return items_[index];
  }

  //! @see operator[](const unsigned)
  PURITY_WEAK constexpr T& at(const unsigned index) noexcept {
    // Not strong purity because items_ is just a pointer!
    return this->operator[](index);
  }

  //! @overload
  PURITY_WEAK constexpr const T& at(const unsigned index) const noexcept {
    // Not strong purity because items_ is just a pointer!
    return this->operator[](index);
  }

  /*! @brief Accessor for the front element
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr T& front() noexcept {
    return items_[0];
  }

  /*! @brief Const accessor for the front element
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr const T& front() const noexcept {
    return items_[0];
  }

  /*! @brief Accessor for the back element
   *
   * @complexity{@math{\Theta(1)}}
   * @note If this is empty, returns front()
   */
  constexpr T& back() noexcept {
    /* No UB in constexpr functions allowed, so we must return something within
     * the array, which is always an initialized value
     */
    if(count_ == 0) {
      return front();
    }

    return items_[count_ - 1];
  }

  /*! @brief Const-accessor for the back element
   *
   * @complexity{@math{\Theta(1)}}
   * @note If this is empty, returns front()
   */
  PURITY_WEAK constexpr const T& back() const noexcept {
    /* NO UB in constexpr functions allowed, so we must return something within
     * the array, which is always an initialized value
     */
    if(count_ == 0) {
      return front();
    }

    return items_[count_ - 1];
  }
//!@}

//!@name Modifiers
//!@{
  /*! @brief Sets the count of contained elements to zero
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr void clear() {
    count_ = 0;
  }

  /*! @brief Copy in elements from another dynamic array
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr void copyIn(const DynamicArray<T, nItems>& other) {
    if(other.size() + size() > nItems) {
      throw "DynamicArray to be copied in has too many elements to fit!";
    }

    for(auto it = other.begin(); it != other.end(); ++it) {
      push_back(*it);
    }
  }

  /*! @brief Insert an element at a particular position
   *
   * @complexity{linear in the number of elements between the inserted position
   * and the end}
   */
  constexpr void insertAt(
    iterator insertPositionIter,
    const T& item
  ) {
    // In case there is no object larger than the one being inserted, add to end
    if(insertPositionIter == end()) {
      push_back(item);
    } else {
      moveElementsRightUntil_(insertPositionIter);

      // Copy in the item
      *insertPositionIter = item;
    }
  }

  /*! @brief Move inserts an element at a specified position
   *
   * @complexity{linear in the number of positions between the insert position
   * and the end}
   */
  constexpr void insertAt(
    iterator insertPositionIter,
    T&& item
  ) {
    // In case there is no object larger than the one being inserted, add to end
    if(insertPositionIter == end()) {
      push_back(item);
    } else {
      moveElementsRightUntil_(insertPositionIter);

      // Move in the item
      *insertPositionIter = std::move(item);
    }
  }

  /*! @brief Removes an item at a particular iterator position
   *
   * @complexity{linear in the number of positions right of the deleted
   * position}
   */
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
//!@}

//!@name Operators
//!@{
  /*! @brief Lexicographical equality comparison
   *
   * @complexity{@math{O(N)}}
   */
  PURITY_WEAK constexpr bool operator == (const DynamicArray& other) const noexcept {
    if(count_ != other.count_) {
      return false;
    }

    for(unsigned i = 0; i < count_; ++i) {
      if(items_[i] != other.items_[i]) {
        return false;
      }
    }

    return true;
  }

  //! Inverts equality comparison
  PURITY_WEAK constexpr bool operator != (const DynamicArray& other) const noexcept {
    return !(*this == other);
  }

  /*! @brief Lexicographical less than comparison
   *
   * @complexity{@math{O(N)}}
   */
  PURITY_WEAK constexpr bool operator < (const DynamicArray& other) const noexcept {
    if(count_ < other.count_) {
      return true;
    }

    for(unsigned i = 0; i < count_; ++i) {
      if(items_[i] < other.items_[i]) {
        return true;
      }

      if(items_[i] > other.items_[i]) {
        return false;
      }
    }

    return false;
  }

  //! Inverts less-than comparison
  PURITY_WEAK constexpr bool operator > (const DynamicArray& other) const noexcept {
    return (other < *this);
  }
//!@}

//!@name Iterators
//!@{
  /**
   * @brief Modifiable data iterator
   */
  class iterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = T&;

    constexpr explicit iterator(
      DynamicArray& instance,
      unsigned&& initPosition
    ) : baseRef_(instance),
        position_(initPosition)
    {}

    constexpr iterator(const iterator& other)
      : baseRef_(other.baseRef_),
        position_(other.position_)
    {}

    constexpr iterator& operator = (const iterator& other) {
      baseRef_ = other.baseRef_;
      position_ = other.position_;

      return *this;
    }

    constexpr iterator& operator ++ () {
      position_ += 1;
      return *this;
    }

    constexpr iterator operator ++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr iterator& operator --() {
      position_ -= 1;
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
      position_ += increment;
      return *this;
    }

    constexpr iterator& operator -= (const int& increment) {
      position_ -= increment;
      return *this;
    }

    PURITY_WEAK constexpr std::ptrdiff_t operator - (const iterator& other) const noexcept {
      return (
        static_cast<std::ptrdiff_t>(position_)
        - static_cast<std::ptrdiff_t>(other.position_)
      );
    }

    PURITY_WEAK constexpr bool operator == (const iterator& other) const noexcept {
      return (
        &baseRef_ == &other.baseRef_
        && position_ == other.position_
      );
    }

    PURITY_WEAK constexpr bool operator != (const iterator& other) const noexcept {
      return !(
        *this == other
      );
    }

    PURITY_WEAK constexpr reference operator * () const noexcept {
      return baseRef_[position_];
    }

  private:
    DynamicArray& baseRef_;
    unsigned position_;
  };

  PURITY_WEAK constexpr iterator begin() noexcept {
    return iterator(*this, 0);
  }

  PURITY_WEAK constexpr iterator end() noexcept {
    return iterator(*this, count_);
  }

  /**
   * @brief Nonmodifiable data iterator
   */
  class const_iterator {
  private:
    const DynamicArray& baseRef_;
    unsigned position_;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;

    constexpr explicit const_iterator(
      const DynamicArray& instance,
      unsigned&& initPosition
    ) : baseRef_(instance),
        position_(initPosition)
    {}

    constexpr const_iterator(const const_iterator& other)
      : baseRef_(other.baseRef_),
        position_(other.position_)
    {}

    constexpr const_iterator& operator = (const const_iterator& other) {
      if(baseRef_ != other.baseRef_) {
        throw "Trying to assign const_iterator to other base DynamicArray!";
      }

      position_ = other.position_;

      return *this;
    }

    constexpr const_iterator& operator ++ () {
      position_ += 1;
      return *this;
    }

    constexpr const_iterator operator ++ (int) {
      const_iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr const_iterator& operator --() {
      position_ -= 1;
      return *this;
    }

    constexpr const_iterator operator -- (int) {
      const_iterator retval = *this;
      --(*this);
      return retval;
    }

    constexpr const_iterator operator + (const int& increment) {
      const_iterator retval = *this;
      retval += increment;
      return retval;
    }

    constexpr const_iterator operator - (const int& increment) {
      const_iterator retval = *this;
      retval -= increment;
      return retval;
    }

    constexpr const_iterator& operator += (const int& increment) {
      position_ += increment;
      return *this;
    }

    constexpr const_iterator& operator -= (const int& increment) {
      position_ -= increment;
      return *this;
    }

    PURITY_WEAK constexpr std::ptrdiff_t operator - (const const_iterator& other) const noexcept {
      return (
        static_cast<std::ptrdiff_t>(position_)
        - static_cast<std::ptrdiff_t>(other.position_)
      );
    }

    PURITY_WEAK constexpr bool operator == (const const_iterator& other) const noexcept {
      return (
        &baseRef_ == &other.baseRef_
        && position_ == other.position_
      );
    }

    PURITY_WEAK constexpr bool operator != (const const_iterator& other) const noexcept {
      return !(
        *this == other
      );
    }

    PURITY_WEAK constexpr reference operator * () const noexcept {
      return baseRef_[position_];
    }
  };

  PURITY_WEAK constexpr const_iterator begin() const noexcept {
    return const_iterator(*this, 0);
  }

  PURITY_WEAK constexpr const_iterator end() const noexcept {
    return const_iterator(*this, count_);
  }
//!@}

//!@name Converting operators
//!@{
  //! Convert the dynamic array to an array
  PURITY_WEAK constexpr operator std::array<T, nItems> () const noexcept {
    return makeArray(std::make_index_sequence<nItems>{});
  }
//!@}

private:
//!@name State
//!@{
  T items_[nItems];
  std::size_t count_ = 0;
//!@}

//!@name Private member functions
//!@{
  template<std::size_t ... Inds>
  std::array<T, nItems> makeArray(std::index_sequence<Inds...> /* inds */) {
    return {{
      items_[Inds]...
    }};
  }

  constexpr void moveElementsRightUntil_(const iterator& position) {
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
//!@}
};

/*!
 * Groups data by equality.
 */
template<
  template<typename, std::size_t> class ArrayType,
  typename T,
  std::size_t N,
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

template<typename T, std::size_t N>
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

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
