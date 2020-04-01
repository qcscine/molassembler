/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Helper class to create strongly typed indices
 */
#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_STRONG_INDEX_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_STRONG_INDEX_H

#include "molassembler/Temple/OperatorSuppliers.h"
#include "boost/functional/hash.hpp"
#include <tuple>

namespace Scine {
namespace temple {

/**
 * @brief Type helper for creating strong index types that are type-level
 *   distinct from their fundamental types
 *
 * @tparam Tag Distinct tag type for each distinct 'type' of tag to avoid
 *   interconvertibility
 * @tparam T Fundamental type of the index
 *
 * Example: Defining a new strong index type
 * @code{.cpp}
 * struct foo_tag;
 * using Foo = StrongIndex<foo_tag, unsigned>;
 * Foo foo(4);
 * @endcode
 *
 * @note This "strong" index type isn't yet as strong as it could be...
 * Implicit base conversion operators aren't great.
 */
template<typename Tag, typename T>
class StrongIndex : public crtp::LexicographicComparable<StrongIndex<Tag, T>> {
public:
  using value_type = T;
  using Hash = boost::hash<StrongIndex<Tag, T>>;

  //! Default constructor does not value initialize
  constexpr StrongIndex() = default;
  //! Explicit value initialization
  constexpr explicit StrongIndex(T v) : v_(v) {}
  constexpr StrongIndex(const StrongIndex& other) = default;
  constexpr StrongIndex(StrongIndex&& other) = default;
  constexpr StrongIndex& operator = (const StrongIndex& other) noexcept = default;
  constexpr StrongIndex& operator = (StrongIndex&& other) noexcept = default;

  //! Assignment from fundamental base type
  constexpr StrongIndex& operator = (const T& v) { v_ = v; return *this; }

  /* Implicit conversion operators */
  constexpr operator T& () noexcept { return v_; }
  constexpr operator const T& () const noexcept { return v_; }

  constexpr auto tie() const { return std::tie(v_); }

private:
  T v_;
};

/**
 * @brief Hash support function for strong index types
 *
 * Allows the use of strong index types in unordered STL containers like so:
 * @code{.cpp}
 * struct foo_tag;
 * using Foo = StrongIndex<foo_tag, unsigned>;
 * using UnorderedSetType = std::unordered_set<Foo, boost::hash<Foo>>;
 * @endcode
 *
 */
template<typename Tag, typename T>
std::size_t hash_value(const StrongIndex<Tag, T>& v) {
  return boost::hash<T>{}(v);
}

} // namespace temple
} // namespace Scine

#endif
