/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides an identity functor.
 */

#ifndef INCLUDE_TEMPLE_FUNCTOR_H
#define INCLUDE_TEMPLE_FUNCTOR_H

#include <algorithm>
#include <stdexcept>
#include "temple/Binding.h"

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
struct At : public Binding<Container> {
  using ContainerBinding = Binding<Container>;
  using ContainerBinding::ContainerBinding;

  template<typename T>
  auto operator() (const T& index) const {
    return ContainerBinding::value.at(index);
  }
};

template<class Container>
auto at(Container&& container) {
  return At<Container&&> {std::forward<Container>(container)};
}

template<class Container>
struct IndexIn : public Binding<Container> {
  using ContainerBinding = Binding<Container>;
  using ContainerBinding::ContainerBinding;

  template<typename T>
  auto operator() (const T& item) const {
    auto findIter = std::find(
      std::begin(ContainerBinding::value),
      std::end(ContainerBinding::value),
      item
    );

    if(findIter == std::end(ContainerBinding::value)) {
      throw std::out_of_range("Item not found in container");
    }

    return findIter - std::begin(ContainerBinding::value);
  }
};

template<class Container>
auto indexIn(Container&& container) {
  return IndexIn<Container&&> {std::forward<Container>(container)};
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
