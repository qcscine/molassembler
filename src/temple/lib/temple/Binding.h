/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides an identity functor.
 */

#ifndef INCLUDE_TEMPLE_BINDING_H
#define INCLUDE_TEMPLE_BINDING_H

#include <utility>

namespace Scine {
namespace temple {

/**
 * @brief Type that will own rvalues, reference lvalues
 *
 * @tparam T Type to conditionally bind
 *
 * In some situations, templated arguments need to be owned by the recipient.
 * For instance, if you want to create a new range from some container, the
 * handling of whether the container being modified is owned by the new range
 * is very important to avoid both extraneous copies and dangling references.
 *
 * We have to differentiate between rvalue and lvalue reference arguments to
 * the adaptor constructor. Any rvalues have to be owned by the range being
 * constructed, while any lvalues should be bound to a (const) reference.
 *
 * This class makes handling that stuff a little easier by providing a single
 * constructor that takes care of all cases.
 *
 * You can also just use this class as a metafunction (see type).
 */
template<typename T>
struct Binding {
  /* Since T is in a type-deducing context, T&& is a 'universal reference',
   * and reference collapsing is applied!
   */
  using type = std::conditional_t<
    std::is_rvalue_reference<T&&>::value || std::is_fundamental<std::decay_t<T>>::value,
    std::decay_t<T>,
    const T&
  >;

  // Now the constructor can apply the same logic unconditionally
  constexpr explicit Binding(T&& t) noexcept : value(t) {}
  virtual ~Binding() = default;

  type value;
};

} // namespace temple
} // namespace Scine

#endif
