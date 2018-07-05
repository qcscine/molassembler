#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_ADD_TO_CONTAINER_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_ADD_TO_CONTAINER_H

#include "ContainerTraits.h"
#include "Functor.h"

/*! @file
 *
 * This file provides the functionality to add an element to a container by
 * calling any of the following members, as appropriate:
 *
 * - insert
 * - emplace
 * - push_back
 * - emplace_back
 */

namespace temple {

/*!
 * Adds an element to a container by calling the container's insert or push_back
 * member functions, if one of either exists.
 */
template<class Container, typename T>
void addToContainer(
  Container& container,
  const T& value
);

/*!
 * Adds an element to a container by calling the container's emplace or
 * emplace_back member functions, if one of either exists.
 */
template<class Container, typename T>
void addToContainer(
  Container& container,
  T&& value
);

template<
  class TargetContainer,
  class SourceContainer,
  class SizeModifierUnary = Identity
> void reserveIfPossible(
  TargetContainer& target,
  const SourceContainer& source,
  SizeModifierUnary&& sourceModifierUnary = Identity()
);

namespace detail {

// SFINAE insertOrPushBack
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

template<class Container, typename T>
std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    T
  >::value && traits::hasEmplaceBack<Container>::value,
  void
> emplaceOrEmplaceBack(
  Container& container,
  const T& value
) {
  container.emplace_back(value);
}

template<class Container, typename T>
std::enable_if_t<
  std::is_same<
    traits::getValueType<Container>,
    T
  >::value && traits::hasEmplace<Container>::value,
  void
> emplaceOrEmplaceBack(
  Container& container,
  const T& value
) {
  container.emplace(value);
}

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
  class SizeModifierUnary = Identity
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
