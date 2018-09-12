#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_TINY_SETS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_TINY_SETS_H

#include <vector>
#include <algorithm>

/*! @file
 *
 * @brief Vector-adapting set-like objects for small collections
 *
 * Vector-based set-like objects for very small collections in order to reduce
 * space overhead and avoid <set>'s trees for better performance.
 */

namespace temple {

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
  using UnderlyingContainer = std::vector<T>;
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
//!@}

//!@name State
//!@{
  UnderlyingContainer data;
//!@}

//!@name Constructors
//!@{
  TinyUnorderedSet() {
    static_assert(
      std::is_same<T, typename std::decay<T>::type>::value,
      "Tiny sets can only contain non-cv qualified types"
    );
  }

  TinyUnorderedSet(std::initializer_list<T> list)
    : data(std::forward<std::initializer_list<T>>(list))
  {
    static_assert(
      std::is_same<T, typename std::decay<T>::type>::value,
      "Tiny sets can only contain non-cv qualified types"
    );
  }
//!@}

//!@name Modification
//!@{
  void insert(T a) {
    data.push_back(std::move(a));
  }

  void emplace(T a) {
    data.push_back(std::move(a));
  }

  template<typename It>
  void erase(It a) {
    auto backIterator = --end();

    if(a == backIterator) {
      data.erase(a);
    } else {
      std::swap(*a, *backIterator);
      data.erase(backIterator);
    }
  }

  template<typename It>
  void insert(It a, const It& b) {
    while(a != b) {
      if(count(*a) == 0) {
        insert(*a);
      }
      ++a;
    }
  }

  void reserve(std::size_t size) {
    data.reserve(size);
  }
//!@}

//!@name Information
//!@{
  unsigned count(T a) const {
    return std::find(
      std::begin(data),
      std::end(data),
      a
    ) != std::end(data);
  }

  std::size_t size() const {
    return data.size();
  }

  bool empty() const {
    return size() == 0;
  }
//!@}

//!@name Iterators
//!@{
  iterator begin() {
    return std::begin(data);
  }

  iterator end() {
    return std::end(data);
  }

  const_iterator begin() const {
    return std::begin(data);
  }

  const_iterator end() const {
    return std::end(data);
  }

  const_iterator cbegin() const {
    return std::cbegin(data);
  }

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
  using UnderlyingContainer = std::vector<T>;
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
//!@}

//!@name State
//!@{
  UnderlyingContainer data;
//!@}

//!@name Constructors
//!@{
  TinySet() {
    static_assert(
      std::is_same<T, typename std::decay<T>::type>::value,
      "Tiny sets can only contain non-cv qualified types"
    );
  }
  TinySet(std::initializer_list<T> list)
    : data(std::forward<std::initializer_list<T>>(list))
  {
    static_assert(
      std::is_same<T, typename std::decay<T>::type>::value,
      "Tiny sets can only contain non-cv qualified types"
    );
  }
//!@}

//!@name Modification
//!@{
  void clear() {
    data.clear();
  }

  template<typename It>
  void erase(It a) {
    data.erase(a);
  }

  void reserve(std::size_t size) {
    data.reserve(size);
  }

  iterator find(const T& a) {
    return std::lower_bound(
      std::begin(data),
      std::end(data),
      a
    );
  }

  void insert(T a) {
    data.insert(
      find(a),
      std::move(a)
    );
  }

  template<typename It>
  void insert(It a, const It& b) {
    while(a != b) {
      if(count(*a) == 0) {
        insert(*a);
      }
      ++a;
    }
  }

  void emplace(T&& a) {
    insert(std::move(a));
  }

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
  unsigned count(T a) const {
    return std::binary_search(
      std::begin(data),
      std::end(data),
      a
    );
  }

  bool empty() const {
    return size() == 0;
  }

  std::size_t size() const {
    return data.size();
  }
//!@}

//!@name Element access
//!@{
  const_reference front() const {
    return data.front();
  }

  const_reference back() const {
    return data.back();
  }

  const_reference at(std::size_t position) const {
    return data.at(position);
  }
//!@}

//!@name Iterators
//!@{
  iterator begin() {
    return std::begin(data);
  }

  iterator end() {
    return std::end(data);
  }

  const_iterator begin() const {
    return std::begin(data);
  }

  const_iterator end() const {
    return std::end(data);
  }

  const_iterator cbegin() const {
    return std::cbegin(data);
  }

  const_iterator cend() const {
    return std::cend(data);
  }
//!@}

//!@name Operators
//!@{
  bool operator == (const TinySet& other) const {
    return data == other.data;
  }

  bool operator != (const TinySet& other) const {
    return !(*this == other);
  }
//!@}
};

} // namespace temple

#endif
