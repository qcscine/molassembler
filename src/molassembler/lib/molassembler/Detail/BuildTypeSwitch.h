/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Defines a compile-time constant boolean indicating the build type.
 */

#ifndef INCLUDE_MOLASSEMBLER_BUILD_TYPE_SWITCH_H
#define INCLUDE_MOLASSEMBLER_BUILD_TYPE_SWITCH_H

namespace molassembler {

#ifdef NDEBUG
constexpr bool buildTypeIsDebug = false;
#else
constexpr bool buildTypeIsDebug = true;
#endif

} // namespace molassembler

#endif
