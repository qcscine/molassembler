/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Container for partially ordered sets with support for gradual ordering
 *   discovery
 */

#ifndef INCLUDE_TEMPLE_POSET_H
#define INCLUDE_TEMPLE_POSET_H

#include "boost/optional.hpp"
#include <vector>
#include <tuple>
#include <bitset>

namespace temple {

/*!
 * @brief Data structure for partially ordered sets with support for gradual
 *   ordering discovery
 */
template<typename T>
class Poset {
public:
  struct Subset {
    using const_iterator = typename std::vector<T>::const_iterator;

    std::vector<T> values;
    bool ordered = false;

    //! Single-value constructor
    Subset(T value) {
      values.push_back(std::move(value));
    }

    //! Vector constructor
    Subset(std::vector<T> passValues) : values(std::move(passValues)) {}

    //! Range constructor
    template<typename Iter>
    Subset(Iter begin, Iter end) : values(begin, end) {}

    const_iterator begin() const {
      return values.begin();
    }

    const_iterator end() const {
      return values.end();
    }
  };

  using const_iterator = typename std::vector<Subset>::const_iterator;

  //! Container constructor
  template<typename Container>
  Poset(const Container& container) {
    setUnorderedValues(
      std::begin(container),
      std::end(container)
    );
  }

  //! Range constructor
  template<typename Iter>
  Poset(Iter&& begin, Iter&& end) {
    setUnorderedValues(
      std::forward<Iter>(begin),
      std::forward<Iter>(end)
    );
  }

  /*!
   * @brief Sets the values in the specified range as the only unordered Subset
   * @tparam Iter The type of the iterator specifying the range
   * @param begin The begin iterator of the range of values
   * @param end The end iterator of the range of values
   */
  template<typename Iter>
  void setUnorderedValues(Iter&& begin, Iter&& end) {
    subsets.clear();
    subsets.emplace_back(std::forward<Iter>(begin), std::forward<Iter>(end));
  }

  /*!
   * @brief Adds order to unordered Subsets
   *
   * @tparam Comparator A callable object with the signature (const T&, const
   *   T&) -> bool
   * @param comparator An instance of Comparator, which returns true if the
   *   first argument is ordered less than the second argument. If two elements
   *   cannot be distinguished with the ordering imposed by @p comparator, then
   *   they should be equal, i.e. !(a < b) && !(b < a).
   */
  template<typename Comparator>
  void orderUnordered(Comparator&& comparator) {
    auto lowerBoundComparator = [&comparator](const Subset& set, const T& t) -> bool {
      return comparator(set.values.front(), t);
    };

    std::vector<Subset> newSubsets;
    newSubsets.reserve(subsets.size());
    for(Subset& subset : subsets) {
      // If the subset is already ordered, move it to the newSubsets
      if(subset.ordered) {
        newSubsets.push_back(std::move(subset));
        continue;
      }

      // Unordered Subsets may split into multiple Subsets
      std::vector<Subset> splat;
      for(auto& value : subset) {
        // Find which Subset this element belongs to using binary search
        auto findIter = std::lower_bound(
          std::begin(splat),
          std::end(splat),
          value,
          lowerBoundComparator
        );

        if(findIter == std::end(splat)) {
          // This element does not belong to any existing Subset, so we create one
          splat.emplace_back(std::move(value));
          continue;
        }

        if(!comparator(value, findIter->values.front())) {
          // This element belongs to the set pointed to by findIter
          findIter->values.push_back(std::move(value));
          continue;
        }

        /* Last possibility: This element is part of a new set that should be
         * inserted before the set pointed to by findIter
         */
        splat.emplace(
          findIter,
          std::move(value)
        );
      }

      // Finalize single-element sets
      for(auto& splatSet : splat) {
        if(splatSet.values.size() == 1) {
          splatSet.ordered = true;
        }
      }

      // Move the split up sets into the new subsets
      std::move(
        std::begin(splat),
        std::end(splat),
        std::back_inserter(newSubsets)
      );
    }

    // Overwrite subsets with the new ones
    subsets = std::move(newSubsets);
  }

  /*!
   * @brief Sets all subsets as ordered
   *
   * Once all comparators are exhausted to find order in yet unordered Subsets,
   * then those values that are still indistinguishable should be considered
   * equal.
   */
  void finalize() {
    for(Subset& subset : subsets) {
      subset.ordered = true;
    }
  }

  const_iterator begin() const {
    return std::begin(subsets);
  }

  const_iterator end() const {
    return std::end(subsets);
  }

  std::string toString() const {
    std::string str = "{";
    for(const Subset& subset : subsets) {
      auto it = std::begin(subset);
      while(true) {
        str += std::to_string(*it);
        if(++it != std::end(subset)) {
          str += ", ";
        } else {
          str += " | ";
          break;
        }
      }
    }
    str += "}";

    return str;
  }

private:
  std::vector<Subset> subsets;
};

} // namespace temple

#endif
