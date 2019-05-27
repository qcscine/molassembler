/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides cyclic shifting frame adaptor to a list-like container
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CYCLIC_FRAME_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CYCLIC_FRAME_ADAPTOR_H

#include "temple/ContainerTraits.h"
#include "temple/constexpr/TupleType.h"
#include <array>
#include <cassert>

namespace temple {
namespace adaptors {
namespace detail {

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
struct CyclicFrameAdaptor {
//!@name Types
//!@{
  // See tricks documentation
  using BoundContainer = std::conditional_t<
    std::is_rvalue_reference<Container&&>::value,
    std::decay_t<Container>,
    const Container&
  >;

  using ContainerValueType = traits::getValueType<Container>;

  using FrameType = typename TupleType::RepeatType<ContainerValueType, frameSize>::type;

  using ContainerIteratorType = decltype(
    std::begin(std::declval<const Container>())
  );
//!@}

//!@name Member variables
//!@{
  BoundContainer container;
//!@}

//!@name Information
//!@{
  std::enable_if_t<
    traits::hasSize<Container>::value,
    std::size_t
  > size() const {
    if(container.size() < frameSize) {
      return 0;
    }

    return container.size();
  }
//!@}

//!@name Special member functions
//!@{
  CyclicFrameAdaptor(Container&& passContainer)
    : container(std::forward<Container>(passContainer)) {}
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
    ) : _begin(begin), _iter(std::move(begin)), _end(std::move(end)) {}

    Iterator& operator ++ () {
      ++_iter;
      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return _iter == other._iter;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    FrameType operator * () const {
      std::array<ContainerValueType, frameSize> frame;

      unsigned firstPhaseCopyCount = std::min(
        frameSize,
        static_cast<unsigned>(std::distance(_iter, _end))
      );

      auto copyIter = std::copy(
        _iter,
        _iter + firstPhaseCopyCount,
        std::begin(frame)
      );

      if(copyIter != std::end(frame)) {
        assert(firstPhaseCopyCount < frameSize);
        const unsigned toCopyCount = frameSize - firstPhaseCopyCount;

        std::copy(
          _begin,
          _begin + toCopyCount,
          copyIter
        );
      }

      return splitArray(frame);
    }

  private:
    ContainerIteratorType _begin, _iter, _end;
  };

  Iterator begin() const {
    // If the container is not of appropriate size, yield an end iterator
    if(container.size() < frameSize) {
      return end();
    }

    return {
      std::begin(container),
      std::end(container)
    };
  }

  Iterator end() const {
    return {
      std::end(container),
      std::end(container)
    };
  }
//!@}
};

} // namespace detail

template<unsigned frameSize, class Container>
auto cyclicFrame(Container&& container) {
  return detail::CyclicFrameAdaptor<Container, frameSize>(
    std::forward<Container>(container)
  );
}

} // namespace adaptors
} // namespace temple

#endif
