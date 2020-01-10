/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief String-specific algorithms
 */

#ifndef INCLUDE_TEMPLE_STRING_ALGORITHMS_H
#define INCLUDE_TEMPLE_STRING_ALGORITHMS_H

#include <string>
#include <sstream>
#include <vector>

namespace Scine {
namespace temple {

inline std::vector<std::string> split(const std::string& s, char delimiter) {
  std::vector<std::string> elements;

  std::stringstream ss;
  ss.str(s);
  std::string item;

  while(std::getline(ss, item, delimiter)) {
    elements.push_back(item);
  }

  return elements;
}

} // namespace temple
} // namespace Scine

#endif
