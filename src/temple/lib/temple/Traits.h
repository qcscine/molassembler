/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides a few basic function and container traits.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_TRAITS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_TRAITS_H

#include <vector>

namespace Scine {
namespace temple {

//! @brief Compile-time reflective trait objects
namespace traits {

// Get the base type a container holds via the begin iterator
namespace detail {
  template<class ContainerType>
  struct getValueTypeImpl {
    using type = typename std::remove_const<
      typename std::remove_reference<
        decltype(
          *std::begin(
            std::declval<ContainerType>()
          )
        )
      >::type
    >::type;
  };

  // Specialization for std::vector<bool>, which returns an awkward proxy object
  template<>
  struct getValueTypeImpl<std::vector<bool>> {
    using type = bool;
  };
} // namespace detail

template<class ContainerType>
using getValueType = typename detail::getValueTypeImpl<ContainerType>::type;

template<class Function, typename ...Args>
using functionReturnType = std::result_of_t<Function(Args...)>;

template<class F>
struct FunctionPointerReturnType;

template<class ReturnType, typename ...Args>
struct FunctionPointerReturnType<ReturnType (*)(Args...)> {
  using type = ReturnType;
};

} // namespace traits
} // namespace temple
} // namespace Scine

#endif
