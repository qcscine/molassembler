/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides cyclic shifting frame adaptor to a list-like container
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CYCLIC_FRAME_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CYCLIC_FRAME_ADAPTOR_H

#include "Molassembler/Temple/ContainerTraits.h"
#include "Molassembler/Temple/Binding.h"
#include "Molassembler/Temple/constexpr/TupleType.h"
#include <cassert>

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Adaptors {
namespace Detail {

template<typename T, std::size_t size, std::size_t ... Inds>
auto splitArrayHelper(const std::array<T, size>& arr, std::index_sequence<Inds ...> /* inds */) {
  return std::make_tuple(
    arr[Inds]...
  );
}

template<typename T, std::size_t size>
auto splitArray(const std::array<T, size>& arr) {
  return splitArrayHelper(arr, std::make_index_sequence<size>());
}

template<class Container, unsigned frameSize>
struct CyclicFrameAdaptor : public Binding<Container> {
//!@name Types
//!@{
  using ContainerBinding = Binding<Container>;
  using ContainerValueType = Traits::getValueType<Container>;

  using FrameType = typename Tuples::RepeatType<ContainerValueType, frameSize>::type;

  using ContainerIteratorType = decltype(
    std::begin(std::declval<const Container>())
  );
//!@}

//!@name Information
//!@{
  std::enable_if_t<
    Traits::hasSize<Container>::value,
    std::size_t
  > size() const {
    if(ContainerBinding::value.size() < frameSize) {
      return 0;
    }

    return ContainerBinding::value.size();
  }
//!@}

//!@name Constructors
//!@{
  using ContainerBinding::ContainerBinding;
//!@}

//!@name Iterators
//!@{
  class Iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = FrameType;
    using difference_type = int;
    using pointer = const FrameType*;
    using reference = const FrameType&;

    Iterator() = default;
    Iterator(
      ContainerIteratorType begin,
      ContainerIteratorType end
    ) : begin_(begin), iter_(std::move(begin)), end_(std::move(end)) {}

    Iterator& operator ++ () {
      ++iter_;
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

    FrameType operator * () const {
      std::array<ContainerValueType, frameSize> frame;

      unsigned firstPhaseCopyCount = std::min(
        frameSize,
        static_cast<unsigned>(std::distance(iter_, end_))
      );

      auto copyIter = std::copy(
        iter_,
        iter_ + firstPhaseCopyCount,
        std::begin(frame)
      );

      if(copyIter != std::end(frame)) {
        assert(firstPhaseCopyCount < frameSize);
        const unsigned toCopyCount = frameSize - firstPhaseCopyCount;

        std::copy(
          begin_,
          begin_ + toCopyCount,
          copyIter
        );
      }

      return splitArray(frame);
    }

  private:
    ContainerIteratorType begin_, iter_, end_;
  };

  Iterator begin() const {
    // If the container is not of appropriate size, yield an end iterator
    if(ContainerBinding::value.size() < frameSize) {
      return end();
    }

    return {
      std::begin(ContainerBinding::value),
      std::end(ContainerBinding::value)
    };
  }

  Iterator end() const {
    return {
      std::end(ContainerBinding::value),
      std::end(ContainerBinding::value)
    };
  }
//!@}
};

} // namespace Detail

template<unsigned frameSize, class Container>
auto cyclicFrame(Container&& container) {
  return Detail::CyclicFrameAdaptor<Container, frameSize>(
    std::forward<Container>(container)
  );
}

} // namespace Adaptors
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
