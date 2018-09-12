#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ORDERED_PAIR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ORDERED_PAIR_H

#include "OperatorSuppliers.h"

#include <tuple>
#include <utility>

/*!@file
 *
 * Contains a class imitating std::pair whose member types are homogeneous and
 * ordered.
 */

namespace temple {

/*!
 * @brief A class that imitates std::pair<T, U>, but whose template arguments
 *   are homogeneous and the stored values ordered (.first < .second).
 *
 * @note We can implement iterators for this really cheaply since the standard
 *   guarantees that successively defined same-aligment types are laid out
 *   successively in memory. Since the pair is homogeneous, we can just use
 *   a T pointer as an iterator.
 */
template<typename T>
struct OrderedPair : crtp::AllOperatorsFromTupleMethod<OrderedPair<T>> {
//!@name Types
//!@{
  using iterator = T*;
  using const_iterator = const T*;
//!@}

//!@name State
//!@{
  // Standard guarantees first and second are laid out sucessively in memory
  T first = T{};
  T second = T{};
//!@}

//!@name Constructors
//!@{
  constexpr OrderedPair() = default;

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
  constexpr auto tuple() const {
    return std::tie(first, second);
  }
//!@}

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
