#ifndef INCLUDE_MOLECULEMANIP_BUILD_TYPE_SWITCH_H
#define INCLUDE_MOLECULEMANIP_BUILD_TYPE_SWITCH_H

namespace MoleculeManip {

#ifdef NDEBUG
constexpr bool buildTypeIsDebug = false;
#else
constexpr bool buildTypeIsDebug = true;
#endif

} // namespace MoleculeManip

#endif
