/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides a few basic function and container traits.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_TRAITS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_TRAITS_H

#include <vector>

namespace Scine {
namespace Molassembler {
namespace Temple {

//! @brief Compile-time reflective trait objects
namespace Traits {

// Get the base type a container holds via the begin iterator
namespace Detail {
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
} // namespace Detail

template<class ContainerType>
using getValueType = typename Detail::getValueTypeImpl<ContainerType>::type;

template<class Function, typename ...Args>
using functionReturnType = std::result_of_t<Function(Args...)>;

template<class F>
struct FunctionPointerReturnType;

template<class ReturnType, typename ...Args>
struct FunctionPointerReturnType<ReturnType (*)(Args...)> {
  using type = ReturnType;
};

} // namespace Traits
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
