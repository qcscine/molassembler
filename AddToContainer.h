#ifndef INCLUDE_TEMPLATE_MAGIC_ADD_TO_CONTAINER_H
#define INCLUDE_TEMPLATE_MAGIC_ADD_TO_CONTAINER_H

#include "ContainerTraits.h"

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

namespace TemplateMagic {

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

} // namespace TemplateMagic

#endif
