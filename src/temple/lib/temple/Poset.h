/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Container for partially ordered sets with support for gradual ordering
 *   discovery
 */

#ifndef INCLUDE_TEMPLE_POSET_H
#define INCLUDE_TEMPLE_POSET_H

#include "boost/logic/tribool.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

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
      return std::begin(values);
    }

    const_iterator end() const {
      return std::end(values);
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

  //! Convert the Poset to a string representation
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

  /**
   * @brief Find an element
   *
   * @param a element to search for
   *
   * @return An iterator pointing to the Subset the element is contained in. If
   *   the element is in no Subset of the Poset, the iterator is the end
   *   iterator.
   */
  const_iterator find(const T& a) const {
    return std::find_if(
      std::begin(subsets),
      std::end(subsets),
      [&a](const Subset& subset) -> bool {
        return std::find(
          std::begin(subset),
          std::end(subset),
          a
        ) != std::end(subset);
      }
    );
  }

  /**
   * @brief Less-than compare two elements
   *
   * @param a The first element
   * @param b The second element
   *
   * @throws std::invalid_argument If either element is not in the poset.
   *
   * @return A tribool in the following state:
   *   - true if the first element is less than the second
   *   - false if the second element is equal to or greater than the second
   *   - indeterminate if ordering is unknown
   */
  boost::tribool compare(const T& a, const T& b) const {
    const auto subsetEnd = std::end(subsets);
    auto aIter = find(a);

    if(aIter == subsetEnd) {
      throw std::invalid_argument("The first argument to compare is not in the poset");
    }

    auto bIter = find(b);

    if(bIter == subsetEnd) {
      throw std::invalid_argument("The second argument to compare is not in the poset");
    }

    /* If both elements are not in the same Subset, then they are ordered by
     * their position in the list of Subsets.
     */
    if(aIter != bIter) {
      return aIter < bIter;
    }

    /* So now we know both outer iterators point to the same subset and that
     * neither of those outer iterators is the end iterator. So we can safely
     * look at the Subset pointed to by the (matching) outer iterators.
     *
     * If the subset is ordered, then the elements are equal. Otherwise, the
     * order is indeterminate.
     */
    if(aIter.first->ordered) {
      return false;
    }

    return boost::indeterminate;
  }

private:
  std::vector<Subset> subsets;
};

} // namespace temple

#endif
