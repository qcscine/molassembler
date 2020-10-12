/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides zip iteration for two containers
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ZIP_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ZIP_ADAPTOR_H

#include "Molassembler/Temple/ContainerTraits.h"
#include "Molassembler/Temple/Binding.h"

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Adaptors {
namespace Detail {

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
    Traits::hasSize<ContainerT>::value && Traits::hasSize<ContainerU>::value,
    std::size_t
  > size() const {
    return std::min(containerT.size(), containerU.size());
  }
//!@}

//!@name Iterators
//!@{
  class Iterator {
  private:
    ContainerTIteratorType left_;
    ContainerUIteratorType right_;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = PairType;
    using difference_type = int;
    using pointer = const PairType*;
    using reference = const PairType&;

    Iterator() = default;
    Iterator(ContainerTIteratorType&& left, ContainerUIteratorType&& right)
      : left_ {left},
        right_ {right}
    {}

    Iterator& operator ++ () {
      ++left_;
      ++right_;

      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return left_ == other.left_ && right_ == other.right_;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {
        *left_,
        *right_
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

} // namespace Detail

template<class ContainerT, class ContainerU>
auto zip(ContainerT&& containerT, ContainerU&& containerU) {
  return Detail::Zipper<ContainerT, ContainerU>(
    std::forward<ContainerT>(containerT),
    std::forward<ContainerU>(containerU)
  );
}

} // namespace Adaptors
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
