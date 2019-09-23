/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Central interface to add elements to many types of containers
 *
 * This file provides the functionality to add an element to a container by
 * calling any of the following members, as appropriate:
 *
 * - insert
 * - emplace
 * - push_back
 * - emplace_back
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ADD_TO_CONTAINER_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ADD_TO_CONTAINER_H

#include "temple/ContainerTraits.h"
#include "temple/Functor.h"

namespace temple {

/*!
 * @brief Adds an lvalue element to a container by calling the container's
 *   insert or push_back member functions, if one of either exists.
 */
template<class Container, typename T>
void addToContainer(
  Container& container,
  const T& value
);

/*!
 * @brief Adds an rvalue to a container by calling the container's emplace or
 *   emplace_back member functions, if one of either exists.
 */
template<class Container, typename T>
void addToContainer(
  Container& container,
  T&& value
);

/*!
 * @brief Rerves the space required in the target container if size can
 *   be determined from the source container
 *
 * If the source container implements a size() const member and the target
 * container implements a reserve() member, the target container reserves space
 * according to the result of the source container size() call modified by the
 * sourceModifierUnary.
 *
 * @tparam TargetContainer the Container type in which space should be reserved
 * @tparam SourceContainer The Container type whose size is used to calculate
 *   required size of the target
 * @tparam SizeModifierUnary A unary callable type that transforms the type
 *   returned by SourceContainer.size(), if applicable
 *
 * @param target The container in which space should be reserved
 * @param source The container whose size should determined for reservation in
 *   the target
 * @param sourceModifierUnary A unary function that transforms the size
 *   determined from the source container. This is for if you expect to generate
 *   more or less elements in the target container.
 */
template<
  class TargetContainer,
  class SourceContainer,
  class SizeModifierUnary = functor::Identity
> void reserveIfPossible(
  TargetContainer& target,
  const SourceContainer& source,
  SizeModifierUnary&& sourceModifierUnary = SizeModifierUnary {}
);

namespace detail {

//! Calls push_back on a container with a passed value if push_back exists
template<class Container, typename T>
std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    T
  >::value && traits::hasPushBack<Container>::value,
  void
> insertOrPushBack(
  Container& container,
  const T& value
) {
  container.push_back(value);
}

//! Calls insert on a container with a passed value if insert exists
template<class Container, typename T>
std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    T
  >::value && traits::hasInsert<Container>::value,
  void
> insertOrPushBack(
  Container& container,
  const T& value
) {
  container.insert(value);
}

//! Calls emplace_back on a container with a passed rvalue if emplace_back exists
template<class Container, typename T>
std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    T
  >::value && traits::hasEmplaceBack<Container>::value,
  void
> emplaceOrEmplaceBack(
  Container& container,
  T&& value
) {
  container.emplace_back(std::forward<T>(value));
}

//! Calls empalce on a container with a passed rvalue if emplace exists
template<class Container, typename T>
std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    T
  >::value && traits::hasEmplace<Container>::value,
  void
> emplaceOrEmplaceBack(
  Container& container,
  T&& value
) {
  container.emplace(std::forward<T>(value));
}

//! Reserves demanded size if all required traits are implemented
template<class TargetContainer, class SourceContainer, class SizeModifierUnary>
std::enable_if_t<
  (
    traits::hasSize<SourceContainer>::value
    && traits::hasReserve<TargetContainer>::value
  ),
  void
> reserve(
  TargetContainer& target,
  const SourceContainer& source,
  SizeModifierUnary&& sizeModifierUnary
) {
  target.reserve(
    sizeModifierUnary(
      source.size()
    )
  );
}

//! Does nothing if either of the required traits are not implemented
template<class TargetContainer, class SourceContainer, class SizeModifierUnary>
std::enable_if_t<
  !(
    traits::hasSize<SourceContainer>::value
    && traits::hasReserve<TargetContainer>::value
  ),
  void
> reserve(
  TargetContainer& /* target */,
  const SourceContainer& /* source */,
  SizeModifierUnary&& /* sizeModifierUnary */
) {
  /* do nothing */
}

} // namespace detail

template<class Container, typename T>
void addToContainer(
  Container& container,
  const T& value
) {
  detail::insertOrPushBack(container, value);
}

template<class Container, typename T>
void addToContainer(
  Container& container,
  T&& value
) {
  detail::emplaceOrEmplaceBack(
    container,
    std::forward<T>(value)
  );
}

template<
  class TargetContainer,
  class SourceContainer,
  class SizeModifierUnary
> void reserveIfPossible(
  TargetContainer& target,
  const SourceContainer& source,
  SizeModifierUnary&& sourceModifierUnary
) {
  return detail::reserve(
    target,
    source,
    std::forward<SizeModifierUnary>(sourceModifierUnary)
  );
}

} // namespace temple

#endif
