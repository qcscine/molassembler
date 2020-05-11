/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief std::array-like class
 *
 * Constexpr fixed-size array to replace std::array in C++14. This class is
 * largely unneeded in C++17 since many std::array members are then marked
 * constexpr.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ARRAY_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ARRAY_H

#include "molassembler/Temple/constexpr/Containers.h"

#include <cstddef>
#include <type_traits>
#include <utility>

namespace Scine {
namespace Molassembler {
namespace Temple {

template<typename T, size_t nItems>
class Array {
public:
//!@name Types
//!@{
  using value_type = T;
//!@}

//!@name Special member functions
//!@{
  /*! @brief Helper copy constructor
   *
   * Constructs using another Array and an index_sequence to directly form the
   * array mem-initializer with a parameter pack expansion
   *
   * @complexity{@math{\Theta(N)}}
   */
  template<size_t ... Inds>
  constexpr Array(
    const Array& other,
    std::index_sequence<Inds...> /* inds */
  ) :items_ {other[Inds]...} {}

  /*!  @brief Copy constructor
   *
   * Constructing from another array is tricky since we're technically not
   * allowed to edit items_ in-class, so we delegate to a helper constructor
   * and directly form the mem-initializer
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr Array(const Array& other)
    : Array(other, std::make_index_sequence<nItems>{}) {}

  /*! @brief Helper move constructor
   *
   * Helper move constructor that directly forms the array mem-initializer with
   * a parameter pack expansion
   *
   * @complexity{@math{\Theta(N)}}
   */
  template<size_t ... Inds>
  constexpr Array(
    Array&& other,
    std::index_sequence<Inds...> /* inds */
  ) :items_ {std::move(other[Inds])...} {}

  /*! @brief Move constructor
   *
   * Moves unfortunately aren't possible for the underlying array as in a
   * pointer swap, but instead must be performed with N moves of the stored
   * types.
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr Array(Array&& other) noexcept
    : Array(std::move(other), std::make_index_sequence<nItems>{}) {}

  /*! @brief Copy assignment operator
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr Array& operator = (const Array& other) {
    for(std::size_t i = 0; i < nItems; ++i) {
      items_[i] = other.items_[i];
    }

    return *this;
  }

  /*! @brief Move assignment operator, moves from other's array
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr Array& operator = (Array&& other) noexcept {
    for(std::size_t i = 0; i < nItems; ++i) {
      items_[i] = std::move(other.items_[i]);
    }

    return *this;
  }

  //! Default destructor
  ~Array() = default;
//!@}

//!@name Converting constructors
//!@{
  //! Delegate std::array ctor, using same trick as copy ctor
  template<size_t ... Inds>
  constexpr Array(
    const std::array<T, nItems>& other,
    std::index_sequence<Inds...> /* inds */
  ) :items_ {other[Inds]...} {}

  /*! @brief Construct from std::array using same trick as copy ctor
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr Array(const std::array<T, nItems>& other)
    : Array(other, std::make_index_sequence<nItems>{})
  {}

  /*! @brief Pack constructor for castable arguments
   *
   * @complexity{@math{\Theta(N)}}
   */
  template<typename ...Args>
  constexpr Array(Args... args)
    : items_ {static_cast<T>(args)...}
  {}
//!@}

//!@name Element acces
//!@{
  /*! @brief Safe access to underlying data
   *
   * Possibility for UB is not "allowed" in constexpr functions, so bounds
   * checking is performed.
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr T& operator[] (const std::size_t index) noexcept {
    return items_[index];
  }

  //! @see operator[](const std::size_t)
  constexpr T& at(const std::size_t index) noexcept {
    return items_[index];
  }

  //! @overload
  constexpr const T& operator[] (const std::size_t index) const noexcept {
    return items_[index];
  }

  //! @overload
  constexpr const T& at(const std::size_t index) const noexcept {
    return items_[index];
  }

  /*! @brief Const accessor for the front element
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr const T& front() const noexcept {
    return items_[0];
  }

  /*! @brief Const Accessor for the back element
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr const T& back() const noexcept {
    return items_[nItems - 1];
  }
//!@}

  constexpr size_t size() const noexcept {
    return nItems;
  }

  /*! @brief Lexicographical equality comparison
   *
   * @complexity{@math{O(N)}}
   */
  constexpr bool operator == (const Array& other) const noexcept {
    for(std::size_t i = 0; i < nItems; ++i) {
      if(items_[i] != other.items_[i]) {
        return false;
      }
    }

    return true;
  }

  //! Inverts equality comparison
  constexpr bool operator != (const Array& other) const noexcept {
    for(std::size_t i = 0; i < nItems; ++i) {
      if(items_[i] != other.items_[i]) {
        return true;
      }
    }

    return false;
  }

  /*! @brief Lexicographical less than comparison
   *
   * @complexity{@math{O(N)}}
   */
  constexpr bool operator < (const Array& other) const noexcept {
    for(std::size_t i = 0; i < nItems; ++i) {
      if(items_[i] < other.items_[i]) {
        return true;
      }

      if(items_[i] > other.items_[i]) {
        return false;
      }
    }

    return false;
  }

  //! Inverts less than comparison
  constexpr bool operator > (const Array& other) const noexcept {
    return (other < *this);
  }

  /**
   * @brief Modifiable data iterator
   */
  class iterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = int;
    using pointer = const T*;
    using reference = T&;

    constexpr iterator(
      Array& instance,
      std::size_t initPosition
    ) noexcept
      : baseRef_(instance),
        position_(initPosition)
    {}

    constexpr iterator() = delete;
    constexpr iterator(const iterator& other) noexcept
      : baseRef_(other.baseRef_),
        position_(other.position_)
    {}
    constexpr iterator(iterator&& other) noexcept = default;
    constexpr iterator& operator = (const iterator& other) noexcept {
      baseRef_ = other.baseRef_;
      position_ = other.position_;

      return *this;
    }
    constexpr iterator& operator = (iterator&& other) noexcept = default;
    ~iterator() = default;

    constexpr iterator& operator ++ () noexcept {
      position_ += 1;
      return *this;
    }

    constexpr iterator operator ++ (int) noexcept {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr iterator& operator --() noexcept {
      position_ -= 1;
      return *this;
    }

    constexpr iterator operator -- (int) noexcept {
      iterator retval = *this;
      --(*this);
      return retval;
    }

    constexpr iterator operator + (const int increment) noexcept {
      iterator retval = *this;
      retval += increment;
      return retval;
    }

    constexpr iterator operator - (const int increment) noexcept {
      iterator retval = *this;
      retval -= increment;
      return retval;
    }

    constexpr iterator& operator += (const int increment) noexcept {
      position_ += increment;
      return *this;
    }

    constexpr iterator& operator -= (const int increment) noexcept {
      position_ -= increment;
      return *this;
    }

    PURITY_WEAK constexpr int operator - (const iterator& other) const noexcept {
      return (
        static_cast<int>(position_)
        - static_cast<int>(other.position_)
      );
    }

    PURITY_WEAK constexpr bool operator == (const iterator& other) const noexcept {
      return (
        &baseRef_ == &other.baseRef_
        && position_ == other.position_
      );
    }

    PURITY_WEAK constexpr bool operator != (const iterator& other) const noexcept {
      return !(
        *this == other
      );
    }

    PURITY_WEAK constexpr reference operator * () const noexcept {
      return baseRef_[position_];
    }

  private:
    Array& baseRef_;
    std::size_t position_;
  };

  PURITY_WEAK constexpr iterator begin() noexcept {
    return iterator(*this, 0);
  }

  PURITY_WEAK constexpr iterator end() noexcept {
    return iterator(*this, nItems);
  }

  /**
   * @brief Nonmodifiable data iterator
   */
  class const_iterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = int;
    using pointer = const T*;
    using reference = const T&;

    constexpr explicit const_iterator(
      const Array& instance,
      std::size_t initPosition
    ) noexcept
      : baseRef_(instance),
        position_(initPosition)
    {}

    constexpr const_iterator() = delete;
    constexpr const_iterator(const const_iterator& other) noexcept
      : baseRef_(other.baseRef_),
        position_(other.position_)
    {}
    constexpr const_iterator(const_iterator&& other) noexcept = default;
    constexpr const_iterator& operator = (const const_iterator& other) {
      if(baseRef_ != other.baseRef_) {
        throw "Trying to assign constIterator to other base Array!";
      }

      position_ = other.position_;

      return *this;
    }
    constexpr const_iterator& operator = (const_iterator&& other) noexcept = default;
    ~const_iterator() = default;

    constexpr const_iterator& operator ++ () noexcept {
      position_ += 1;
      return *this;
    }

    constexpr const_iterator operator ++ (int) noexcept {
      const_iterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr const_iterator& operator --() noexcept {
      position_ -= 1;
      return *this;
    }

    constexpr const_iterator operator -- (int) noexcept {
      const_iterator retval = *this;
      --(*this);
      return retval;
    }

    constexpr const_iterator operator + (const int increment) noexcept {
      const_iterator retval = *this;
      retval += increment;
      return retval;
    }

    constexpr const_iterator operator - (const int increment) noexcept {
      const_iterator retval = *this;
      retval -= increment;
      return retval;
    }

    constexpr const_iterator& operator += (const int increment) noexcept {
      position_ += increment;
      return *this;
    }

    constexpr const_iterator& operator -= (const int increment) noexcept {
      position_ -= increment;
      return *this;
    }

    constexpr int operator - (const const_iterator& other) const noexcept {
      return (
        static_cast<int>(position_)
        - static_cast<int>(other.position_)
      );
    }

    constexpr bool operator == (const const_iterator& other) const noexcept {
      return (
        &baseRef_ == &other.baseRef_
        && position_ == other.position_
      );
    }

    constexpr bool operator != (const const_iterator& other) const noexcept {
      return !(
        *this == other
      );
    }

    constexpr reference operator * () const noexcept {
      return baseRef_[position_];
    }

  private:
    const Array& baseRef_;
    std::size_t position_;
  };

  constexpr const_iterator begin() const noexcept {
    return const_iterator(*this, 0);
  }

  constexpr const_iterator end() const noexcept {
    return const_iterator(*this, nItems);
  }

  //! Implicit conversion operator to a std::array
  constexpr operator std::array<T, nItems> () const {
    return makeArray_(std::make_index_sequence<nItems>{});
  }

  //! Explicit conversion to a std::array
  constexpr std::array<T, nItems> getArray() const {
    return makeArray_(std::make_index_sequence<nItems>{});
  }

private:
  T items_[nItems];

  template<size_t ... Inds>
  std::array<T, nItems> makeArray_(std::index_sequence<Inds...> /* inds */) {
    return {{
      items_[Inds]...
    }};
  }
};

template<typename T, typename... Tail>
constexpr auto makeArray(T head, Tail... tail) -> Array<T, 1 + sizeof...(Tail)> {
  // This way, makeArray must be called with at least one argument
  return { head, tail ... };
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
