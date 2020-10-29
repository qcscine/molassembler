/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides functors for transformations
 */

#ifndef INCLUDE_TEMPLE_FUNCTOR_H
#define INCLUDE_TEMPLE_FUNCTOR_H

#include <algorithm>
#include <stdexcept>
#include "Molassembler/Temple/Binding.h"

namespace Scine {
namespace Molassembler {
namespace Temple {

//! @brief Functors for transformations
namespace Functor {

/**
 * @brief Metafunction required for default arguments
 */
struct Identity {
  //! @brief Returns its arguments (perfect forwarding)
  template<typename U>
  constexpr auto operator()(U&& v) const noexcept {
    return std::forward<U>(v);
  }
};

/*! @brief Metafunction calling at on a bound container
 *
 * @note Rvalues and lvalues are acceptable. Rvalues will be owned by the
 * functor, lvalues by reference.
 */
template<class Container>
struct At : public Binding<Container> {
  //! Type of bound container
  using ContainerBinding = Binding<Container>;
  // Bring constructor into scope
  using ContainerBinding::ContainerBinding;

  //! @brief Return value in container at specified position
  template<typename T>
  auto operator() (const T& index) const {
    return ContainerBinding::value.at(index);
  }
};

/**
 * @brief Make functor calling at on arguments
 * @param container container to bind
 * @return functor invokable with index type
 */
template<class Container>
auto at(Container&& container) {
  return At<Container&&> {std::forward<Container>(container)};
}

//! @brief Metafunction determining the index of an item in a bound container
template<class Container>
struct IndexIn : public Binding<Container> {
  //! Type of bound container
  using ContainerBinding = Binding<Container>;
  // Bring constructor into scope
  using ContainerBinding::ContainerBinding;

  /**
   * @brief Determine index of an item in the container
   * @param item sought item
   * @return index of item in the container
   * @throws std::out_of_range if the item is not found
   */
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

/**
 * @brief Make functor determining the index of an item in a bound container
 * @param container container to bind
 * @return functor invokable with index type
 */
template<class Container>
auto indexIn(Container&& container) {
  return IndexIn<Container&&> {std::forward<Container>(container)};
}

//! @brief Metafunction calling std::get on an argument with a bound index
template<unsigned i>
struct Get {
  //! @brief Call std::get<i> on the argument
  template<typename T>
  constexpr auto operator() (const T& t) const {
    return std::get<i>(t);
  }
};

/**
 * @brief Make functor calling std::get with particular index on an argument
 *
 * @tparam i Index of get to use
 *
 * @return Functor invokable with everything you can call std::get on
 */
template<unsigned i>
constexpr auto get() { return Get<i> {}; }

//! @brief Calls std::get<0> on any argument it is invoked with
constexpr Get<0> first;
//! @brief Calls std::get<1> on any argument it is invoked with
constexpr Get<1> second;

template<typename UnaryF, typename UnaryG>
struct Currier {
  Currier(UnaryF&& outer, UnaryG&& inner) : f(outer), g(inner) {}

  template<typename T>
  auto operator() (T&& t) -> decltype(auto) {
    return f(g(t));
  }

  UnaryF f;
  UnaryG g;
};

template<typename UnaryF, typename UnaryG>
auto curry(UnaryF&& outer, UnaryG&& inner) {
  return Currier<UnaryF, UnaryG>(std::forward<UnaryF>(outer), std::forward<UnaryG>(inner));
}

} // namespace Functor
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
