#ifndef INCLUDE_TEMPLATE_MAGIC_TINY_SET_H
#define INCLUDE_TEMPLATE_MAGIC_TINY_SET_H

#include <vector>
#include <algorithm>

namespace temple {

// Break-even with an unordered_set is at around N = 80 for numeric types
template<typename T>
struct TinyUnorderedSet {
  using type = std::vector<T>;

  type data;

  void insert(T a) {
    data.push_back(a);
  }

  template<typename It>
  void insert(It a, It b) {
    while(a != b) {
      if(count(*a) == 0) {
        insert(*a);
      }
      ++a;
    }
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
};

// Break-even with set is at around N = 200 for numeric types
template<typename T>
struct TinySet {
  using type = std::vector<T>;

  type data;

  void insert(T a) {
    data.insert(
      std::lower_bound(
        data.begin(),
        data.end(),
        a
      ),
      a
    );
  }

  template<typename It>
  void insert(It a, It b) {
    while(a != b) {
      if(count(*a) == 0) {
        insert(*a);
      }
      ++a;
    }
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
};

} // namespace temple

#endif
