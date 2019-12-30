/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides zip iteration for two containers
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ZIP_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ZIP_ADAPTOR_H

#include "temple/ContainerTraits.h"
#include "temple/Binding.h"

namespace temple {

namespace adaptors {

namespace detail {

template<class ContainerT, class ContainerU>
struct Zipper {
//!@name Types
//!@{
  // See tricks documentation
  using BoundContainerT = typename Binding<ContainerT>::type;
  using BoundContainerU = typename Binding<ContainerU>::type;

  using T = decltype(
    *std::begin(
      std::declval<const ContainerT>()
    )
  );
  using U = decltype(
    *std::begin(
      std::declval<const ContainerU>()
    )
  );

  using PairType = std::pair<const T&, const U&>;

  using ContainerTIteratorType = decltype(std::begin(std::declval<const ContainerT>()));
  using ContainerUIteratorType = decltype(std::begin(std::declval<const ContainerU>()));
//!@}

//!@name Reference in order to own possible rvalue containers / ranges
//!@{
  BoundContainerT containerT;
  BoundContainerU containerU;
//!@}

//!@name Constructors
//!@{
  Zipper(
    ContainerT&& t,
    ContainerU&& u
  ) : containerT(std::forward<ContainerT>(t)),
      containerU(std::forward<ContainerU>(u))
  {}
//!@}

//!@name Information
//!@{
  std::enable_if_t<
    (
      traits::hasSize<ContainerT>::value
      && traits::hasSize<ContainerU>::value
    ),
    std::size_t
  > size() const {
    return std::min(
      containerT.size(),
      containerU.size()
    );
  }
//!@}

//!@name Iterators
//!@{
  class Iterator {
  private:
    ContainerTIteratorType _left;
    ContainerUIteratorType _right;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = PairType;
    using difference_type = int;
    using pointer = const PairType*;
    using reference = const PairType&;

    Iterator() = default;
    Iterator(ContainerTIteratorType&& left, ContainerUIteratorType&& right)
      : _left {left},
        _right {right}
    {}

    Iterator& operator ++ () {
      ++_left;
      ++_right;

      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return _left == other._left && _right == other._right;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {
        *_left,
        *_right
      };
    }
  };

  Iterator begin() const {
    return {
      std::begin(containerT),
      std::begin(containerU)
    };
  }

  Iterator end() const {
    auto tSize = std::distance(std::begin(containerT), std::end(containerT));
    auto uSize = std::distance(std::begin(containerU), std::end(containerU));

    if(tSize == uSize) {
      return {
        std::end(containerT),
        std::end(containerU)
      };
    }

    if(tSize < uSize) {
      auto uIter = std::begin(containerU);
      std::advance(uIter, tSize);

      return {
        std::end(containerT),
        std::move(uIter)
      };
    }

    auto tIter = std::begin(containerT);
    std::advance(tIter, uSize);

    return {
      std::move(tIter),
      std::end(containerU)
    };
  }
//!@}
};

} // namespace detail

template<class ContainerT, class ContainerU>
auto zip(ContainerT&& containerT, ContainerU&& containerU) {
  return detail::Zipper<ContainerT, ContainerU>(
    std::forward<ContainerT>(containerT),
    std::forward<ContainerU>(containerU)
  );
}

} // namespace adaptors

} // namespace temple

#endif
