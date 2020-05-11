/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Class for flat maps between strong indices
 */
#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_STRONG_INDEX_MAP_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_STRONG_INDEX_MAP_H

#include <vector>
#include <algorithm>
#include <stdexcept>

namespace Scine {
namespace Temple {

template<typename T, typename U>
class StrongIndexFlatMap {
public:
  using key_type = T;
  using value_type = U;

  using iterator = typename std::vector<U>::iterator;
  using const_iterator = typename std::vector<U>::const_iterator;

  StrongIndexFlatMap() = default;
  explicit StrongIndexFlatMap(std::vector<U> strong) : map_(std::move(strong)) {}
  explicit StrongIndexFlatMap(const std::vector<typename U::value_type>& weak)
    : map_(std::begin(weak), std::end(weak)) {}

  T indexOf(const U u) const {
    auto findIter = std::find(
      std::begin(map_),
      std::end(map_),
      u
    );

    if(findIter == std::end(map_)) {
      throw std::out_of_range("Item not found");
    }

    return T(findIter - std::begin(map_));
  }

  void resize(unsigned newSize) {
    map_.resize(newSize);
  }

  unsigned size() const {
    return map_.size();
  }

  U& at(const T t) {
    return map_.at(t);
  }

  const U& at(const T t) const {
    return map_.at(t);
  }

  bool empty() const {
    return map_.empty();
  }

  void pushIsometric() {
    map_.emplace_back(map_.size());
  }

  void clear() {
    map_.clear();
  }

  StrongIndexFlatMap<U, T> invert() const {
    const unsigned N = map_.size();
    std::vector<T> inverted(N);
    for(unsigned i = 0; i < N; ++i) {
      inverted.at(map_.at(i)) = T(i);
    }
    return StrongIndexFlatMap<U, T>(inverted);
  }

  iterator begin() { return std::begin(map_); }
  iterator end() { return std::end(map_); }
  const_iterator begin() const { return std::begin(map_); }
  const_iterator end() const { return std::end(map_); }

  bool operator == (const StrongIndexFlatMap& other) const {
    return map_ == other.map_;
  }

  bool operator != (const StrongIndexFlatMap& other) const {
    return map_ != other.map_;
  }

private:
  std::vector<U> map_;
};

} // namespace Temple
} // namespace Scine

#endif
