/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Vector-adapting set-like objects for small collections
 *
 * Vector-based set-like objects for very small collections in order to reduce
 * space overhead and avoid `std::set`'s trees for better performance.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_TINY_SETS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_TINY_SETS_H

#include <vector>
#include <algorithm>

namespace Scine {
namespace Temple {

/*!
 * @brief An adapter class for std::vector that acts like an unordered set
 *   and stores its data in the underlying vector in an unordered fashion.
 *
 * Why? It's a lot simpler / smaller than unordered_set and seems to be faster
 * for very small amounts of data. Break-even with an unordered_set for finding
 * a value is at around N = 80 for numeric types, although I certainly do not
 * trust the benchmark much.
 *
 * @tparam The set value type
 */
template<typename T>
struct TinyUnorderedSet {
//!@name Types
//!@{
  //! The type of the underlying container with linear memory
  using UnderlyingContainer = std::vector<T>;
  //! What types this container stores
  using value_type = T;
  //! Reference to a value_type
  using reference = T&;
  //! Const reference to a value type
  using const_reference = const T&;
  //! Iterator type
  using iterator = typename std::vector<T>::iterator;
  //! Const iterator type
  using const_iterator = typename std::vector<T>::const_iterator;
//!@}

//!@name State
//!@{
  //! The owned underlying linear memory container
  UnderlyingContainer data;
//!@}

//!@name Constructors
//!@{
  //! Default constructor
  TinyUnorderedSet() {
    static_assert(
      std::is_same<T, typename std::decay<T>::type>::value,
      "Tiny sets can only contain non-cv qualified types"
    );
  }

  //! Initializer-list constructor
  TinyUnorderedSet(std::initializer_list<T> list)
    : data(std::forward<std::initializer_list<T>>(list))
  {
    static_assert(
      std::is_same<T, typename std::decay<T>::type>::value,
      "Tiny sets can only contain non-cv qualified types"
    );
  }
//!@}

//!@name Static functions
//!@{
  /**
   * @brief Inserts an element only if the element is not already in the set
   *
   * @param set Set to insert the element into
   * @param a Element to insert
   *
   * @return Whether the element was already in the set
   */
  static bool checked_insert(std::vector<T>& set, T a) {
    if(!contains(set, a)) {
      set.push_back(std::move(a));
      return false;
    }

    return true;
  }

  /**
   * @brief Checked emplace of an element into a set
   *
   * In-place constructs a T, then checks if it already exists in the set. If
   * it does, it is removed again.
   *
   * @tparam Args Parameter pack to forward to T's constructor
   * @param set Set to insert the element into
   * @param args Arguments to forward to T's constructor
   */
  template<typename ... Args>
  static void checked_emplace(std::vector<T>& set, Args&& ... args) {
    set.emplace_back(std::forward<Args>(args)...);

    if(
      std::find(
        std::begin(set),
        std::end(set) - 1,
        set.back()
      ) != std::end(set) - 1
    ) {
      set.pop_back();
    }
  }

  /**
   * @brief Linear search for a value
   *
   * @param set Set to search
   * @param value Value to look for
   *
   * @return Whether the value was found
   */
  static bool linear_search(std::vector<T>& set, const T& value) {
    return std::find(
      std::begin(set),
      std::end(set),
      value
    ) != std::end(set);
  }
//!@}

//!@name Modification
//!@{
  //! Empties the set
  void clear() {
    data.clear();
  }

  /*!
   * @brief Inserts an element into the set
   * @warning This does not check if the element already exists.
   */
  void insert(T a) {
    data.push_back(std::move(a));
  }

  /*!
   * @brief Emplaces an element into the set
   * @note You can use this just like @p std::vector::emplace_back
   * @warning This does not check if the element already exists
   */
  template<typename ... Args>
  void emplace(Args&& ... args) {
    data.emplace_back(std::forward<Args>(args)...);
  }

  //! Erases a position in the set
  template<typename It>
  void erase(It a) {
    auto backIterator = --end();

    if(a == backIterator) {
      data.erase(a);
    } else {
      std::iter_swap(a, backIterator);
      data.erase(backIterator);
    }
  }

  //! Inserts all elements contained in a range
  template<typename It>
  void insert(It a, const It& b) {
    while(a != b) {
      if(count(*a) == 0) {
        insert(*a);
      }
      ++a;
    }
  }

  //! Reserves space in memory
  void reserve(std::size_t size) {
    data.reserve(size);
  }
//!@}

//!@name Information
//!@{
  //! Counts the number of occurrences of an element
  unsigned count(T a) const {
    return std::find(
      std::begin(data),
      std::end(data),
      a
    ) != std::end(data);
  }

  //! Returns the number of elements in the set
  std::size_t size() const {
    return data.size();
  }

  //! Returns whether the set is empty
  bool empty() const {
    return size() == 0;
  }
//!@}

//!@name Iterators
//!@{
  //! Yields a begin iterator
  iterator begin() {
    return std::begin(data);
  }

  //! Yields an end iterator
  iterator end() {
    return std::end(data);
  }

  //! Yields a begin const iterator
  const_iterator begin() const {
    return std::begin(data);
  }

  //! Yields an end const iterator
  const_iterator end() const {
    return std::end(data);
  }

  //! Yields a begin const iterator
  const_iterator cbegin() const {
    return std::cbegin(data);
  }

  //! Yields an end const iterator
  const_iterator cend() const {
    return std::cend(data);
  }
//!@}
};

/*!
 * @brief An adapter class for std::vector that acts like an ordered set
 *   and stores its data in the underlying vector in order.
 *
 * Why? It's a lot simpler / smaller than set and seems to be faster
 * for very small amounts of data. Break-even with an unordered_set for finding
 * a value is at around N = 200 for numeric types, although I certainly do not
 * trust the benchmark much.
 *
 * @tparam The set value type
 */
template<typename T>
struct TinySet {
//!@name Types
//!@{
  //! The type of the underlying container with linear memory
  using UnderlyingContainer = std::vector<T>;
  //! What types this container stores
  using value_type = T;
  //! Reference to a value_type
  using reference = T&;
  //! Const reference to a value type
  using const_reference = const T&;
  //! Iterator type
  using iterator = typename std::vector<T>::iterator;
  //! Const iterator type
  using const_iterator = typename std::vector<T>::const_iterator;
//!@}

//!@name State
//!@{
  //! The owned underlying linear memory container
  UnderlyingContainer data;
//!@}

//!@name Constructors
//!@{
  //! Default constructor
  TinySet() {
    static_assert(
      std::is_same<T, typename std::decay<T>::type>::value,
      "Tiny sets can only contain non-cv qualified types"
    );
  }
  //! Initializer-list constructor
  TinySet(std::initializer_list<T> list)
    : data(std::forward<std::initializer_list<T>>(list))
  {
    static_assert(
      std::is_same<T, typename std::decay<T>::type>::value,
      "Tiny sets can only contain non-cv qualified types"
    );
  }
//!@}

//!@name Static functions
//!@{
  static iterator find(std::vector<T>& set, const T& a) {
    return std::lower_bound(
      std::begin(set),
      std::end(set),
      a
    );
  }

  static const_iterator find(const std::vector<T>& set, const T& a) {
    return std::lower_bound(
      std::begin(set),
      std::end(set),
      a
    );
  }

  static void checked_insert(std::vector<T>& set, T a) {
    iterator findIter = find(set, a);

    /* Condition lifted directly from std::binary_search. If the element equals
     * a, then the set already contains the element. In that case, do nothing.
     */
    if(findIter != std::end(set) && !(a < *findIter)) {
      return;
    }

    set.insert(std::move(findIter), std::move(a));
  }

  /**
   * @brief Emplaces an element in an ordered vector if the element is not
   *   contained
   *
   * @note This does inplace-construct the element at the back of the vector,
   *   then searches for an equal element in the ordered part of the vector. If
   *   the element is found, the newly-constructed element is popped back off.
   *   Otherwise, it is rotated to its position. This implies that this
   *   function may trigger reallocation of the vector even though the element
   *   is already part of the vector.
   *
   * @param set An ordered vector to add the element to
   * @tparam Args parameter pack to forward to T's constructor
   * @param args Arguments to forward to T's constructor
   */
  template<typename ... Args>
  static void checked_emplace(std::vector<T>& set, Args&& ... args) {
    // In-place construct the element at the back of the set
    set.emplace_back(std::forward<Args>(args)...);

    iterator lastElementIter = std::end(set) - 1;

    // Search the range before it for an identical element
    iterator findIter = std::lower_bound(
      std::begin(set),
      lastElementIter,
      set.back()
    );

    if(findIter != lastElementIter && !(set.back() < *findIter)) {
      // Remove the element again
      set.pop_back();
    } else {
      // Rotate in the new element
      std::rotate(
        findIter,
        lastElementIter,
        std::end(set)
      );
    }
  }

  /**
   * @brief Binary search for a value in an ordered set
   *
   * @param set Set to search
   * @param value Value to look for
   *
   * @return Whether the value was found
   */
  static bool binary_search(const std::vector<T>& set, const T& value) {
    return std::binary_search(
      std::begin(set),
      std::end(set),
      value
    );
  }
//!@}

//!@name Modification
//!@{
  //! Empties the set
  void clear() {
    data.clear();
  }

  //! Erases an element marked by a position in the set
  template<typename It>
  void erase(It a) {
    data.erase(a);
  }

  //! Reserves space in memory
  void reserve(std::size_t size) {
    data.reserve(size);
  }

  //! Finds an element in the set
  iterator find(const T& a) {
    return std::lower_bound(
      std::begin(data),
      std::end(data),
      a
    );
  }

  /*!
   * @brief Inserts an element into the set
   * @warning Does not check if the element already exists
   */
  void insert(T a) {
    data.insert(
      find(a),
      std::move(a)
    );
  }

  /*!
   * @brief Inserts all elements contained in a range
   * @warning Does not check if the element already exists
   */
  template<typename It>
  void insert(It a, const It& b) {
    while(a != b) {
      if(count(*a) == 0) {
        insert(*a);
      }
      ++a;
    }
  }

  /*!
   * @brief Insert an element
   * @note This does not actually in-place construct. This is a copy and move.
   * @warning Does not check if the element already exists
   */
  void emplace(T a) {
    insert(std::move(a));
  }

  /*!
   * @brief Assign from an initializer list
   * @warning Repeated elements are retained in insertion
   */
  TinySet& operator = (std::initializer_list<T> init) {
    data.clear();
    for(const auto& value : init) {
      insert(value);
    }

    return *this;
  }
//!@}

//!@name Information
//!@{
  //! Count the number of occurrences of an element
  unsigned count(T a) const {
    return std::binary_search(
      std::begin(data),
      std::end(data),
      a
    );
  }

  //! Returns whether the set is empty
  bool empty() const {
    return size() == 0;
  }

  //! Returns the number of contained elements
  std::size_t size() const {
    return data.size();
  }
//!@}

//!@name Element access
//!@{
  //! Yields a const reference to the first element in the set
  const_reference front() const {
    return data.front();
  }

  //! Yields a const reference to the last element in the set
  const_reference back() const {
    return data.back();
  }

  //! Yields a const reference to an element in the set
  const_reference at(std::size_t position) const {
    return data.at(position);
  }
//!@}

//!@name Iterators
//!@{
  //! Yields a begin iterator
  iterator begin() {
    return std::begin(data);
  }

  //! Yields an end iterator
  iterator end() {
    return std::end(data);
  }

  //! Yields a begin const iterator
  const_iterator begin() const {
    return std::begin(data);
  }

  //! Yields an end const iterator
  const_iterator end() const {
    return std::end(data);
  }

  //! Yields a begin const iterator
  const_iterator cbegin() const {
    return std::cbegin(data);
  }

  //! Yields an end const iterator
  const_iterator cend() const {
    return std::cend(data);
  }
//!@}

//!@name Operators
//!@{
  //! Lexicographically compare two sets
  bool operator == (const TinySet& other) const {
    return data == other.data;
  }

  //! Invertes @p operator==
  bool operator != (const TinySet& other) const {
    return !(*this == other);
  }
//!@}
};

} // namespace Temple
} // namespace Scine

#endif
