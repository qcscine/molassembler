#ifndef INCLUDE_MOLASSEMBLER_VERSION_H
#define INCLUDE_MOLASSEMBLER_VERSION_H

#include <string>

/*! @file
 *
 * Contains the library versioning scheme information
 */

namespace molassembler {

namespace Version {

constexpr unsigned short major = 0;
constexpr unsigned short minor = 1;
constexpr unsigned short fix = 0;

inline std::string majorMinor() {
  return std::to_string(major) + "." + std::to_string(minor);
}

inline std::string fullVersion() {
  return (
    std::to_string(major)
    + "." + std::to_string(minor)
    + "." + std::to_string(fix)
  );
}

} // namespace Version

} // namespace molassembler

#endif
