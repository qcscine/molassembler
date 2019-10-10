/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides an identity functor.
 */

#ifndef INCLUDE_TEMPLE_IDENTITY_H
#define INCLUDE_TEMPLE_IDENTITY_H

#include <utility>

namespace temple {

namespace functor {

/**
 * @brief Metafunction required for default arguments
 */
struct Identity {
  template<typename U>
  constexpr auto operator()(U&& v) const noexcept {
    return std::forward<U>(v);
  }
};

/*! @brief Metafunction calling at on any arguments
 *
 * @note Rvalues and lvalues are acceptable. Rvalues will be owned by the
 * functor, lvalues by reference.
 */
template<class Container>
struct At {
  using BoundContainer = std::conditional_t<
    std::is_rvalue_reference<Container&&>::value,
    std::decay_t<Container>,
    const Container&
  >;

  BoundContainer bound;

  At(Container&& container) : bound(container) {}

  template<typename T>
  auto operator() (const T& index) const {
    return bound.at(index);
  }
};

template<class Container>
auto at(Container&& container) {
  return At<Container&&> {std::forward<Container>(container)};
}

template<unsigned i>
struct Get {
  template<typename T>
  constexpr auto operator() (const T& t) const {
    return std::get<i>(t);
  }
};

template<unsigned i>
constexpr auto get() { return Get<i> {}; }

constexpr Get<0> first;
constexpr Get<1> second;

} // namespace functor

} // namespace temple

#endif
