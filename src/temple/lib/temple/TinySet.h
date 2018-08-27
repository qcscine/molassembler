#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_TINY_SET_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_TINY_SET_H

#include <vector>
#include <algorithm>

/*! @file
 *
 * Vector-based set-like objects for very small collections in order to reduce
 * space overhead and better performance.
 */

namespace temple {

// Break-even with an unordered_set is at around N = 80 for numeric types
template<typename T>
struct TinyUnorderedSet {
  using type = std::vector<T>;

  type data;

  TinyUnorderedSet() = default;
  TinyUnorderedSet(std::initializer_list<T> list) : data {std::forward<std::initializer_list<T>>(list)} {}

  void insert(T a) {
    data.push_back(std::move(a));
  }

  void emplace(T a) {
    data.push_back(std::move(a));
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

  unsigned count(T a) const {
    return std::find(
      data.begin(),
      data.end(),
      a
    ) != data.end();
  }

  unsigned size() const {
    return data.size();
  }

  bool empty() const {
    return size() == 0;
  }

  typename type::iterator begin() {
    return data.begin();
  }

  typename type::iterator end() {
    return data.end();
  }

  typename type::const_iterator begin() const {
    return data.begin();
  }

  typename type::const_iterator end() const {
    return data.end();
  }

  typename type::const_iterator cbegin() const {
    return data.cbegin();
  }

  typename type::const_iterator cend() const {
    return data.cend();
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
};

// Break-even with set is at around N = 200 for numeric types
template<typename T>
struct TinySet {
  using type = std::vector<T>;

  type data;

  TinySet() = default;
  TinySet(std::initializer_list<T> list) : data {std::forward<std::initializer_list<T>>(list)} {}

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

  typename type::iterator find(const T& a) {
    return std::lower_bound(
      data.begin(),
      data.end(),
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

  unsigned size() const {
    return data.size();
  }

  unsigned count(T a) const {
    return std::binary_search(
      data.begin(),
      data.end(),
      a
    );
  }

  T front() const {
    return data.front();
  }

  T back() const {
    return data.back();
  }

  T at(std::size_t position) const {
    return data.at(position);
  }

  bool empty() const {
    return size() == 0;
  }

  typename type::iterator begin() {
    return data.begin();
  }

  typename type::iterator end() {
    return data.end();
  }

  typename type::const_iterator begin() const {
    return data.begin();
  }

  typename type::const_iterator end() const {
    return data.end();
  }

  typename type::const_iterator cbegin() const {
    return data.cbegin();
  }

  typename type::const_iterator cend() const {
    return data.cend();
  }

  bool operator == (const TinySet& other) const {
    return data == other.data;
  }

  bool operator != (const TinySet& other) const {
    return !(*this == other);
  }
};

} // namespace temple

#endif
