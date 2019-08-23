/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief An ordered homogeneous pair-like class
 *
 * Contains a class imitating std::pair whose member types are homogeneous and
 * ordered.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ORDERED_PAIR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ORDERED_PAIR_H

#include "OperatorSuppliers.h"

#include <tuple>
#include <utility>

namespace temple {

/*!
 * @brief A class that imitates std::pair<T, U>, but whose template arguments
 *   are homogeneous and the stored values ordered (.first < .second).
 *
 * @note We can implement iterators for this really cheaply since the standard
 *   guarantees that successively defined same-aligment types are laid out
 *   successively in memory. Since the pair is homogeneous, we can just use
 *   a T pointer as an iterator. These do not, however, fulfill all the
 *   requirements for STL algorithms such as having typedefs for the results of
 *   various operations.
 */
template<typename T>
struct OrderedPair : crtp::LexicographicComparable<OrderedPair<T>> {
//!@name Types
//!@{
  //! Type of stored elements
  using value_type = T;
  //! Iterator type
  using iterator = T*;
  //! Const iterator type
  using const_iterator = const T*;
//!@}

//!@name State
//!@{
  //! First element of the pair
  T first;
  //! Second element of the pair
  T second;
//!@}

//!@name Constructors
//!@{
  //! Default constructor
  constexpr OrderedPair() = default;

  //! Reordering pair initializer
  constexpr OrderedPair(T a, T b) : first {std::move(a)}, second {std::move(b)} {
    if(b < a) {
      std::swap(first, second);
    }
  }
//!@}

//!@name Element access
//!@{
  constexpr T front() const {
    return first;
  }

  constexpr T& front() {
    return first;
  }

  constexpr T back() const {
    return second;
  }

  constexpr T& back() {
    return second;
  }
//!@}

//!@name Iterators
//!@{
  iterator begin() {
    return &first;
  }

  iterator end() {
    return std::next(&second);
  }

  const_iterator begin() const {
    return &first;
  }

  const_iterator end() const {
    return std::next(&second);
  }

  const_iterator cbegin() const {
    return &first;
  }

  const_iterator cend() const {
    return std::next(&second);
  }
//!@}

//!@name Operators
//!@{
  //! Yields the result of std::tie(first, second)
  constexpr auto tie() const {
    return std::tie(first, second);
  }
//!@}

  //! Map the pair into an unordered std::pair
  template<typename UnaryFunction>
  auto map(UnaryFunction&& mapFunction) const {
    return std::make_pair(
      mapFunction(first),
      mapFunction(second)
    );
  }
};

} // namespace temple

#endif
