#ifndef INCLUDE_MOLECULEMANIP_BUILD_TYPE_SWITCH_H
#define INCLUDE_MOLECULEMANIP_BUILD_TYPE_SWITCH_H

/*! @file
 * Defines a compile-time constant boolean indicating the build type.
 */

namespace molassembler {

#ifdef NDEBUG
constexpr bool buildTypeIsDebug = false;
#else
constexpr bool buildTypeIsDebug = true;
#endif

} // namespace molassembler

#endif
