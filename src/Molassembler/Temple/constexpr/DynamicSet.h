/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief BTree-based std::set-like container (but max size is space allocated)
 *
 * A constexpr fixed-maximum-size managed set so that the type signature does
 * not change upon element insertion and deletion. STL parallel is std::set,
 * with the difference that here the maximum number of elements must be known at
 * compile time.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_DYNAMIC_SET_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_DYNAMIC_SET_H

#include "Molassembler/Temple/constexpr/BTree.h"

namespace Scine {
namespace Molassembler {
namespace Temple {

/**
 * @brief Tree-based set
 *
 * @tparam T Value type of the set
 * @tparam nItems Maximum number of items to store in the set
 * @tparam LessThanPredicate
 */
template<
  typename T,
  size_t nItems,
  class LessThanPredicate = std::less<T>,
  class EqualityPredicate = std::equal_to<T>
> class DynamicSet {
private:
  using TreeType = BTree<T, 3, nItems, LessThanPredicate, EqualityPredicate>;

  TreeType tree_;

public:
  //! Dynamic set
  constexpr DynamicSet() = default;

  /*! @brief Constructor from existing ordered data
   *
   * Warning: These constructors expect ordered arrays!
   */
  template<
    template<typename, size_t> class ArrayType,
    size_t size
  > constexpr DynamicSet(const ArrayType<T, size>& items) {
    for(const auto& item : items) {
      tree_.insert(item);
    }
  }

  /*! @brief Check if the set contains an element.
   *
   * @complexity{@math{\Theta(N \log N)}}
   */
  PURITY_WEAK constexpr bool contains(const T& item) const {
    return tree_.contains(item);
  }

  /*! @brief Insertion an element into the set.
   *
   * @complexity{@math{\Theta(N \log N)}}
   */
  constexpr void insert(const T& item) {
    tree_.insert(item);
  }

  PURITY_WEAK constexpr Optional<const T&> getOption(const T& item) const {
    return tree_.getOption(item);
  }

  //! @brief Remove all elements of the set
  constexpr void clear() {
    tree_.clear();
  }

  using const_iterator = typename TreeType::const_iterator;

  PURITY_WEAK constexpr const_iterator begin() const {
    return tree_.begin();
  }

  PURITY_WEAK constexpr const_iterator end() const {
    return tree_.end();
  }

  PURITY_WEAK constexpr size_t size() const {
    return tree_.size();
  }

  PURITY_WEAK constexpr bool operator == (const DynamicSet& other) const {
    return tree_ == other.tree_;
  }

  PURITY_WEAK constexpr bool operator != (const DynamicSet& other) const {
    return !(
      tree_ == other.tree_
    );
  }

  PURITY_WEAK constexpr bool operator < (const DynamicSet& other) const {
    return tree_ < other.tree_;
  }

  PURITY_WEAK constexpr bool operator > (const DynamicSet& other) const {
    return other.tree_ < tree_;
  }
};

//! Helper function to create a DynamicSet specifying only the maximum size
template<
  size_t nItems,
  typename T,
  template<typename, size_t> class ArrayType,
  class LessThanPredicate = std::less<T>,
  class EqualityPredicate = std::equal_to<T>
> constexpr DynamicSet<T, nItems, LessThanPredicate> makeDynamicSet(
  const ArrayType<T, nItems>& array
) {
  return DynamicSet<T, nItems, LessThanPredicate, EqualityPredicate>(array);
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
