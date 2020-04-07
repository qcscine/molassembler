/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Defines a compile-time constant boolean indicating the build type.
 */

#ifndef INCLUDE_MOLASSEMBLER_BUILD_TYPE_SWITCH_H
#define INCLUDE_MOLASSEMBLER_BUILD_TYPE_SWITCH_H

namespace Scine {

namespace molassembler {

#ifdef NDEBUG
constexpr bool buildTypeIsDebug = false;
#else
constexpr bool buildTypeIsDebug = true;
#endif

} // namespace molassembler

} // namespace Scine

#endif
