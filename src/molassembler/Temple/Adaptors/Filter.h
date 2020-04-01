/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides filtering to ranges with arbitrary predicates
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_FILTER_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_FILTER_ADAPTOR_H

#include "temple/ContainerTraits.h"
#include "temple/Binding.h"
#include "temple/constexpr/TupleType.h"
#include <cassert>

namespace Scine {
namespace temple {
namespace adaptors {
namespace detail {

template<class Container, typename UnaryPredicate>
struct FilterAdaptor {
//!@name Types
//!@{
  // See tricks documentation
  using BoundContainer = typename Binding<Container>::type;

  using ContainerValueType = traits::getValueType<Container>;

  using ContainerIteratorType = decltype(
    std::begin(std::declval<const Container>())
  );
//!@}

//!@name Member variables
//!@{
  BoundContainer container;
  UnaryPredicate predicate;
//!@}

//!@name Special member functions
//!@{
  FilterAdaptor(Container&& passContainer, UnaryPredicate&& passPredicate)
    : container(std::forward<Container>(passContainer)), predicate(passPredicate) {}
//!@}

//!@name Iterators
//!@{
  class Iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ContainerValueType;
    using difference_type = int;
    using pointer = const ContainerValueType*;
    using reference = const ContainerValueType&;

    Iterator() = default;
    Iterator(
      ContainerIteratorType begin,
      ContainerIteratorType end,
      const FilterAdaptor& parent
    ) : iter_(std::move(begin)), end_(std::move(end)), adaptor_(parent) {
      if(iter_ != end_ && !adaptor_.get().predicate(*iter_)) {
        do {
          ++iter_;
        } while(iter_ != end_ && !adaptor_.get().predicate(*iter_));
      }
    }

    Iterator& operator ++ () {
      do {
        ++iter_;
      } while(iter_ != end_ && !adaptor_.get().predicate(*iter_));
      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return iter_ == other.iter_;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    ContainerValueType operator * () const {
      return *iter_;
    }

  private:
    ContainerIteratorType iter_, end_;
    std::reference_wrapper<const FilterAdaptor> adaptor_;
  };

  Iterator begin() const {
    return {
      std::begin(container),
      std::end(container),
      *this
    };
  }

  Iterator end() const {
    return {
      std::end(container),
      std::end(container),
      *this
    };
  }
//!@}
};

} // namespace detail

/**
 * @brief Keeps elements satifying a predicate (i.e. "keep filter")
 *
 * @param container The container to adapt
 * @param predicate The predicate specifying which elements to keep
 *
 * @return An adapted range object
 */
template<class Container, typename UnaryPredicate>
auto filter(Container&& container, UnaryPredicate&& predicate) {
  return detail::FilterAdaptor<Container, UnaryPredicate>(
    std::forward<Container>(container),
    std::forward<UnaryPredicate>(predicate)
  );
}

} // namespace adaptors
} // namespace temple
} // namespace Scine

#endif
